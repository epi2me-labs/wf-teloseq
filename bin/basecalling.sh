#!/usr/bin/env bash
set -e

usage() {
    echo "Usage: $0 -m <model> -i <pod5_directory> -o <output_directory>"
    echo "Valid models: 'hac', 'sup'"
    exit 1
}

while getopts ":m:i:o:" opt; do
    case ${opt} in
        m) model=$OPTARG ;;
        i) input=$OPTARG ;;
        o) output=$OPTARG ;;
        *) usage ;;
    esac
done

if [[ "$model" != "hac" && "$model" != "sup" ]]; then
    usage
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

dorado basecaller \
	--barcode-sequences "${tmp_dir}/teloseq_adapters.fasta" \
	--barcode-arrangement "${tmp_dir}/teloseq.toml" \
	--kit-name telo-seq --no-trim \
	"${model}" "${input}" \
	| dorado demux --no-classify --output-dir "${output}"

echo "Preparing output structure..."
pushd "${output}"
for bc in {01..12}; do
    mkdir "barcode$bc"
    mv ./*barcode"$bc".bam "barcode$bc"
done
mkdir unclassified
mv ./*unclassified.bam unclassified
popd
