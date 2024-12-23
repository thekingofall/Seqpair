import gzip
import re
import argparse
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Function to create a regex pattern for the reverse complement with ambiguous bases
def reverse_complement_ambiguous_to_regex(seq):
    complement = {
        'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C',
        'M': '[KT]', 'S': '[CS]',
        'K': '[MT]', 'R': '[YA]', 'Y': '[RG]', 'B': '[VGT]', 'V': '[BCT]',
        'D': '[HAT]', 'H': '[DBT]', 'N': '.'
    }
    return ''.join(complement.get(base, base) for base in reversed(seq))

# Function to search for patterns with up to 2 mismatches
def search_with_mismatches(pattern, text, max_mismatches=2):
    for m in re.finditer(pattern, text, re.IGNORECASE):
        if sum(1 for a, b in zip(m.group(), text[m.start():m.end()]) if a != b) <= max_mismatches:
            return True
    return False

# Argument parsing
parser = argparse.ArgumentParser(description="Search for specific sequences in paired-end FASTQ files using regex.")
parser.add_argument("basename", help="The base name for the FASTQ files (without _R1/_R2 or file extension)")
args = parser.parse_args()

# Define the sequences and create their regex patterns
seq1 = 'GAGGACCTGAAMAASGTGT'
seq2 = 'ATATCCAGAACCCTGACCC'
seq1_rc_regex = reverse_complement_ambiguous_to_regex(seq1)
seq2_rc_regex = reverse_complement_ambiguous_to_regex(seq2)

# Initialize counters
read_count = 0
match_count = 0

# Process the FASTQ files
with gzip.open(f"{args.basename}_R1.fastq.gz", 'rt') as handle1, \
     gzip.open(f"{args.basename}_R2.fastq.gz", 'rt') as handle2:
    for (title1, seq1, qual1), (title2, seq2, qual2) in zip(FastqGeneralIterator(handle1), FastqGeneralIterator(handle2)):
        read_count += 1
        if read_count > 1000:  # Stop after 1000 reads
            break
        seq1_found = search_with_mismatches(seq1_rc_regex, seq1)
        seq2_found = search_with_mismatches(seq2_rc_regex, seq2)
        if seq1_found and seq2_found:
            match_count += 1

# Calculate and print the percentage
percentage = (match_count / read_count) * 100 if read_count else 0
print(f"Total reads processed: {read_count}")
print(f"Reads with both sequences: {match_count} ({percentage:.2f}%)")