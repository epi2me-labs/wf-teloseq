#!/usr/bin/env bash
set -euo pipefail


usage() {
  cat <<USAGE
Usage: $0 [-m <model>] [-d <path to dorado>] -i <pod5_directory> -o <output_directory>
Arguments:
  -i  Input POD5 directory
  -o  Output directory
Options:
  -m  Dorado model (hac|sup)
  -d  Dorado executable (defaults to 'dorado' in PATH)
USAGE
  exit 1
}


dorado="dorado"
model="hac"
input=""
output=""

while getopts ":hm:i:o:d:" opt; do
    case ${opt} in
        h) usage ;;
        m) model=$OPTARG ;;
        i) input=$OPTARG ;;
        o) output=$OPTARG ;;
        d) dorado=$OPTARG ;;
        *) echo "Error: Unexpected option '$opt'."
	   echo
	   usage ;;
    esac
done

[[ -n "${input}" ]] || { echo "Error: -i <pod5_directory> is required." >&2; echo; usage; }
[[ -n "${output}" ]] || { echo "Error: -o <output_directory> is required." >&2; echo; usage; }

if [[ "$model" != "hac" && "$model" != "sup" ]]; then
    echo "Error: Invalid model '$model'. Valid models: 'hac', 'sup'" >&2
    echo
    usage
fi

if [[ ! -x "$dorado" ]]; then
    if ! command -v "$dorado" >/dev/null 2>&1 ; then
	    echo "Error: Cannot execute '$dorado'. Please consider using the -d option to specify its location." >&2
	    echo
        usage
    fi
fi

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

cat << EOF > "$tmp_dir/teloseq_adapters.fasta"
>TA01
CACAAAGACACCGACAACTTTCTT
>TA02
AAGGTTACACAAACCCTGGACAAG
>TA03
AAGGATTCATTCCCACGGTAACAC
>TA04
GAGAGGACAAAGGTTTCAACGCTT
>TA05
TCCGATTCTGCTTCTTTCTACCTG
>TA06
AGAACGACTTCCATACTCGTGTGA
>TA07
CGTCAACTGACAGTGGTTCGTACT
>TA08
CCAAACCCAACAACCTAGATAGGC
>TA09
CCAGTAGAAGTCCGACAACGTCAT
>TA10
GGAGTTCGTCCAGAGAAGTACACG
>TA11
CTTTCGTTGTTGACTCGACGGTAG
>TA12
CATCTGGAACGTGGTACACCTGTA
EOF

cat << EOF > "$tmp_dir/teloseq.toml"
[arrangement]
name = "telo-seq"
kit = "EXP-TLO001"

mask1_front = "ATTGCTAAGGTTAA"
mask1_rear = "CCCTAACC"

# Barcode sequences
barcode1_pattern = "TA%02i"
first_index = 1
last_index = 12
EOF

echo "Using dorado executable: $dorado"
# Get dorado version (e.g. "1.3.0+6ea400189")
version_str=$("$dorado" --version 2>&1)
echo "Dorado version string: ${version_str:-unknown}"
version_num=$(printf '%s\n' "$version_str" | grep -oE '[0-9]+\.[0-9]+\.[0-9]+' | head -n1)
echo "Dorado version: ${version_num:-unknown}"

# dorado >= 1.2.0 cannot read demux input from stdin, and has a different output format
demux_needs_file=false
if [ -n "$version_num" ]; then
    if [ "$(printf '1.2.0\n%s\n' "$version_num" | sort -V | head -n1)" = "1.2.0" ]; then
        demux_needs_file=true
    fi
fi

# dorado demux --no-classify doesn't work in v1.2.0
if [ "$version_num" = "1.2.0" ]; then
    echo "Error: dorado v1.2.0 has a known issue with 'demux --no-classify'. Please upgrade dorado to v1.3.0 or later." >&2
    exit 1
fi

echo "Creating output directory: ${output}"
mkdir -p "${output}"

basecaller_cmd=(
    "$dorado" basecaller
    --barcode-sequences "${tmp_dir}/teloseq_adapters.fasta"
    --barcode-arrangement "${tmp_dir}/teloseq.toml"
    --kit-name telo-seq
    --no-trim
    "${model}" "${input}"
)
echo "Running basecalling and demultiplexing..."

if [ "$demux_needs_file" = false ]; then
    echo "Using dorado version $version_num: demux can read from stream"
    "${basecaller_cmd[@]}" \
        | "$dorado" demux --no-classify --output-dir "${output}"

    echo "Preparing output structure..."
    pushd "${output}"
    for bc in {01..12}; do
        mkdir "barcode$bc"
        mv ./*barcode"$bc".bam "barcode$bc"
    done
    mkdir unclassified
    mv ./*unclassified.bam unclassified
    popd
else
    echo "Using dorado version $version_num: demux requires intermediate file"
    tmp_bam="${output}/mixed.tmp.bam"
    "${basecaller_cmd[@]}" > "${tmp_bam}"

    "$dorado" demux \
        --no-classify \
        --output-dir "${output}" \
        "${tmp_bam}"
    rm -f "${tmp_bam}"

    # output for later dorado versions is MinKNOW like
    echo "Preparing output structure..."
    pushd "${output}"
    mv ./*/*/bam_pass/* .
    for bc in {01..12}; do
        mv telo-seq_barcode"$bc" "barcode$bc"
    done
    # rm anything not barcodeXX or unclassified
    find . -maxdepth 1 ! -name 'barcode*' ! -name 'unclassified' ! -name '.' -exec rm -rf {} +
fi
