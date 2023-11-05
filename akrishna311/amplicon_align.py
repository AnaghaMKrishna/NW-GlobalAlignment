#!/usr/bin/env python3

import magnumopus
import argparse

def reverse_complement(string:str) -> str:
    """
    Generate the complement of the string using complement_dict

    Args:
        string : string to be complemented

    Returns:
        The complemented string
    """
    complement_dict = {
    "A":"T",
    "G":"C",
    "T":"A",
    "C":"G"
}
    return "".join(complement_dict.get(base) for base in string[::-1])

parser = argparse.ArgumentParser(prog="amplicon_align.py", 
                                 description="Perform in-silico PCR on two assemblies and align the amplicons using Needleman-Wunsch global alignment algorithm")

arg_assembly_file1 = parser.add_argument(
    "-1", "--assembly-file1",
    help="Path to the first assembly file",
    dest="assembly1",
    type=str,
    required=True
)

arg_assembly_file2 = parser.add_argument(
    "-2", "--assembly-file2",
    help="Path to the second assembly file",
    dest="assembly2",
    type=str,
    required=True
)

arg_primer_file = parser.add_argument(
    "-p", "--primer-file",
    help="Path to the primer file",
    dest="primers",
    type=str,
    required=True
)

arg_ampl_size = parser.add_argument(
    "-m", "--max-amplicon-size",
    help="maximum amplicon size for isPCR",
    dest="max_amplicon_size",
    type=int,
    required=True
)

arg_match = parser.add_argument(
    "--match",
    help="match score to use in alignment",
    dest="match",
    type=int,
    required=True
)

arg_mismatch = parser.add_argument(
    "--mismatch",
    help="mismatch penalty to use in alignment",
    dest="mismatch",
    type=int,
    required=True
)

arg_gap = parser.add_argument(
    "--gap",
    help="gap penalty to use in alignment",
    dest="gap",
    type=int,
    required=True
)

# parser.print_help()
args = parser.parse_args()

assembly_file1 = args.assembly1
assembly_file2 = args.assembly2
primer_file = args.primers
max_amplicon_size = args.max_amplicon_size
match = args.match
mismatch = args.mismatch
gap = args.gap

amplicon1 = magnumopus.ispcr(primer_file, assembly_file1, max_amplicon_size)
# print(amplicon1)
amplicon2 = magnumopus.ispcr(primer_file, assembly_file2, max_amplicon_size)
# print(amplicon2)
amplicon1 = amplicon1.lstrip(amplicon1[:amplicon1.index("\n")]).strip()
amplicon2 = amplicon2.lstrip(amplicon2[:amplicon2.index("\n")]).strip()

# print(f"{amplicon1}\n\n\n{amplicon2}")


aln1, score1 = magnumopus.needleman_wunsch(amplicon1, amplicon2, match, mismatch, gap)
aln2, score2 = magnumopus.needleman_wunsch(reverse_complement(amplicon1), amplicon2, match, mismatch, gap)
aln3, score3 = magnumopus.needleman_wunsch(amplicon1, reverse_complement(amplicon2), match, mismatch, gap)
aln4, score4 = magnumopus.needleman_wunsch(reverse_complement(amplicon1), reverse_complement(amplicon2), match, mismatch, gap)
max_score = [[score1, aln1], [score2, aln2], [score3, aln3], [score4, aln4]]
opt_aln = max_score.index(max(max_score))
print("\n".join(max_score[opt_aln][1]))
print(max_score[opt_aln][0])

