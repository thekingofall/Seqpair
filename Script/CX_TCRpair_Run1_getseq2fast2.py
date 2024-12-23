import gzip
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import regex
import sys
from tqdm import tqdm  # Import tqdm for the progress bar

# Check if the correct number of arguments are provided
if len(sys.argv) != 4:
    print("Usage: python script.py <input_fastq_basename> <sequence_to_search> <max_mismatches>")
    sys.exit(1)

input_basename = sys.argv[1]
search_sequence = sys.argv[2].upper()
max_mismatches = int(sys.argv[3])

# Function to create a regex pattern allowing for ambiguous bases and mismatches
def generate_pattern(sequence, max_mismatches):
    sequence = (sequence.replace('S', '[CG]')
                        .replace('M', '[AC]'))
    # Allow up to max_mismatches
    return regex.compile(f"({sequence}){{e<={max_mismatches}}}")

# Function to open the fastq.gz files
def open_fastq(path):
    return gzip.open(path, 'rt')

# Generate the patterns for search
forward_pattern = generate_pattern(search_sequence, max_mismatches)
reverse_pattern = generate_pattern(search_sequence[::-1], max_mismatches)
complement_pattern = generate_pattern(str(Seq(search_sequence).complement()), max_mismatches)
rev_complement_pattern = generate_pattern(str(Seq(search_sequence).reverse_complement()), max_mismatches)

# Counters for the occurrences
# ... (the rest of your imports and initial code)

# Counters for the occurrences
forward_count = reverse_count = complement_count = rev_complement_count = non_match_count = 0
total_count = 0  # Total records counter

# Read the fastq files and search for the patterns
with open_fastq(f"{input_basename}_1.fq.gz") as fastq_r1, open_fastq(f"{input_basename}_2.fq.gz") as fastq_r2:
    for (title_r1, sequence_r1, quality_r1), (title_r2, sequence_r2, quality_r2) in zip(FastqGeneralIterator(fastq_r1), FastqGeneralIterator(fastq_r2)):
        total_count += 1  # Increment total records counter
        match_found = False
        
        # Search in R1
        if forward_pattern.search(sequence_r1):
            forward_count += 1
            match_found = True
        if reverse_pattern.search(sequence_r1):
            reverse_count += 1
            match_found = True
        if complement_pattern.search(sequence_r1):
            complement_count += 1
            match_found = True
        if rev_complement_pattern.search(sequence_r1):
            rev_complement_count += 1
            match_found = True

        # Search in R2
        if forward_pattern.search(sequence_r2):
            forward_count += 1
            match_found = True
        if reverse_pattern.search(sequence_r2):
            reverse_count += 1
            match_found = True
        if complement_pattern.search(sequence_r2):
            complement_count += 1
            match_found = True
        if rev_complement_pattern.search(sequence_r2):
            rev_complement_count += 1
            match_found = True

        # Increment non-match count if no match was found in either read
        if not match_found:
            non_match_count += 1

        # Stop after processing 100 lines (50 records)
        if total_count >= 10000:
            break

# Calculate percentages
forward_percentage = (forward_count / total_count * 100) if total_count > 0 else 0
reverse_percentage = (reverse_count / total_count * 100) if total_count > 0 else 0
complement_percentage = (complement_count / total_count * 100) if total_count > 0 else 0
rev_complement_percentage = (rev_complement_count / total_count * 100) if total_count > 0 else 0
non_match_percentage = (non_match_count / total_count * 100) if total_count > 0 else 0

# Print the results
print(f"Total reads processed: {total_count}")
print(f"Forward sequence count: {forward_count} ({forward_percentage:.2f}%)")
print(f"Reverse sequence count: {reverse_count} ({reverse_percentage:.2f}%)")
print(f"Complement sequence count: {complement_count} ({complement_percentage:.2f}%)")
print(f"Reverse complement sequence count: {rev_complement_count} ({rev_complement_percentage:.2f}%)")
print(f"Non-matching reads count: {non_match_count} ({non_match_percentage:.2f}%)")