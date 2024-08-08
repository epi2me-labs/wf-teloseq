"""Extract reference from FASTA file."""

import sys

import pysam

from .util import get_named_logger, wf_parser  # noqa: ABS101


def reverse_complement(seq):
    """Get the reverse complement of a DNA sequence."""
    bases = {
        "A": "T",
        "C": "G",
        "G": "C",
        "T": "A",
    }
    revcomp = ""
    try:
        for base in seq[::-1]:
            revcomp += bases[base]
    except KeyError as e:
        (unexpected_char,) = e.args
        raise ValueError(f"Found unexpected character '{unexpected_char}' in sequence.")
    return revcomp


def extract_telomere(
    seq,
    telomere_seq,
    restriction_motif,
    extra_spacing,
    reverse=False,
):
    """Extract telomere + sub-telomere on forward or reverse strand of a sequence."""
    # find first occurrence of restriction site
    seq = seq.upper()
    if reverse:
        seq = reverse_complement(seq)
    cut_pos = seq.find(restriction_motif.upper())
    if cut_pos != -1:
        # we found a cut site; let's make sure that there is a telomere upstream
        if telomere_seq.upper() in seq[:cut_pos]:
            end_pos = cut_pos + extra_spacing
            return seq[:end_pos]
    # no telomere found
    return None


def main(args):
    """Run the entry point."""
    logger = get_named_logger("extractRefs")

    with pysam.FastxFile(args.fasta_file) as f:
        for entry in f:
            # we first check the forward and then the reverse strand (`from_end`) to
            # make sure that there is a telomere and restriction site
            for reverse, suffix in ((False, "_p"), (True, "_q")):
                telomere = extract_telomere(
                    entry.sequence,
                    telomere_seq=args.telomere_seq,
                    restriction_motif=args.restriction_motif,
                    extra_spacing=args.extra_spacing,
                    reverse=reverse,
                )
                if telomere is None:
                    orientation = "reverse" if reverse else "forward"
                    logger.info(f"No {orientation} telomere found for '{entry.name}'.")
                    continue
                # check if the ref sequence ID ends in "_p" or "_q" (designating the
                # short or long chromosome arm, respecitvely; if not, we assume that the
                # short arm comes first)
                corrected_header = (
                    entry.name
                    if entry.name[-2:].lower() in ("_p", "_q")
                    else f"{entry.name}{suffix}"
                )
                sys.stdout.write(f">{corrected_header}\n{telomere}\n")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("extract_reference")
    parser.add_argument(
        "fasta_file",
        help="Input reference FASTA file",
    )
    parser.add_argument(
        "telomere_seq",
        help="Telomere sequence to search for",
    )
    parser.add_argument(
        "restriction_motif",
        help="Motif of restriction enzyme",
    )
    parser.add_argument(
        "extra_spacing",
        type=int,
        help="Number of extra basepairs to include in the extracted sequence",
    )
    return parser
