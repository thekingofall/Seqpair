import gzip
import re
import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def check_patterns(seq):
    patterns = [
        re.compile(r'GGGTCAGGGTTCTGGATAT', re.IGNORECASE),
        re.compile(r'AGCTGCTATGCACGACTG', re.IGNORECASE),
        re.compile(r'GAGGACCTGAA', re.IGNORECASE),
        re.compile(r'GTGT', re.IGNORECASE)
    ]
    return all(pattern.search(seq) for pattern in patterns)

def extract_umis_and_remove_pattern(seq):
    pattern = r'GAGGACCTGA[AC][AC]AA[CG]GTGT([ACGT]{8})AGCTGCTATGCACGACTG([ACGT]{8})GGGTCAGGGTTCTGGATAT'
    match = re.search(pattern, seq, flags=re.IGNORECASE)
    if not match:
        return None, seq
    umi1 = match.group(1)
    umi2 = match.group(2)
    umi = f"UMI:{umi1}{umi2}"
    
    # Remove the matched pattern from the sequence
    new_seq = seq[:match.start()] + seq[match.end():]
    
    return umi, new_seq

def extract_umis_and_remove_pattern(seq):
    pattern = r'GAGGACCTGA[AC][AC]AA[CG]GTGT([ACGT]{8})AGCTGCTATGCACGACTG([ACGT]{8})GGGTCAGGGTTCTGGATAT'
    match = re.search(pattern, seq, flags=re.IGNORECASE)
    if not match:
        return None, seq
    umi1 = match.group(1)
    umi2 = match.group(2)
    umi = f"UMI:{umi1}{umi2}"
    
    # Remove the matched pattern from the sequence
    new_seq = seq[:match.start()] + seq[match.end():]
    
    return umi, new_seq
def process_fastq_pair(r1_file, r2_file, out_file_prefix, barcode_file, max_reads):
    out_r1_file = out_file_prefix + '_R1.fastq.gz'
    out_r2_file = out_file_prefix + '_R2.fastq.gz'
    
    with gzip.open(r1_file, 'rt') as f1, gzip.open(r2_file, 'rt') as f2, \
         gzip.open(out_r1_file, 'wt') as out1, gzip.open(out_r2_file, 'wt') as out2, \
         open(barcode_file, 'w') as barcode_out:
        
        r1_records = SeqIO.parse(f1, 'fastq')
        r2_records = SeqIO.parse(f2, 'fastq')
        
        barcode_out.write("read_id\tumi\tread_origin\n")
        
        matched_reads = 0
        total_reads = 0
        debug_count = 0

        for r1, r2 in zip(r1_records, r2_records):
            total_reads += 1
            
            r1_seq = str(r1.seq)
            r2_seq = str(r2.seq)

            r1_match = check_patterns(r1_seq)
            r2_match = check_patterns(r2_seq)

            if debug_count < 5:
                print(f"\nDebug info for read pair {debug_count + 1}:")
                print(f"R1 sequence: {r1_seq}")
                print(f"R2 sequence: {r2_seq}")
                print(f"R1 patterns match: {r1_match}")
                print(f"R2 patterns match: {r2_match}")
                debug_count += 1

            umi = None
            if r1_match:
                umi, remaining_seq = extract_umis_and_remove_pattern(r1_seq)
                if umi:
                    matched_reads += 1
                    if matched_reads <= 5:
                        print(f"Matched R1: {r1.description}")
                    new_seq = Seq(remaining_seq)
                    new_qual = r1.letter_annotations["phred_quality"][:len(remaining_seq)]
                    r1 = SeqRecord(new_seq, id=r1.id, name=r1.name, description=r1.description)
                    r1.letter_annotations["phred_quality"] = new_qual
            elif r2_match:
                umi, remaining_seq = extract_umis_and_remove_pattern(r2_seq)
                if umi:
                    matched_reads += 1
                    if matched_reads <= 5:
                        print(f"Matched R2: {r2.description}")
                    new_seq = Seq(remaining_seq)
                    new_qual = r2.letter_annotations["phred_quality"][:len(remaining_seq)]
                    r2 = SeqRecord(new_seq, id=r2.id, name=r2.name, description=r2.description)
                    r2.letter_annotations["phred_quality"] = new_qual
            
            if umi:
                r1.description += f" {umi}"
                r2.description += f" {umi}"
                SeqIO.write(r1, out1, 'fastq')
                SeqIO.write(r2, out2, 'fastq')
                barcode_out.write(f"{r1.id}\t{umi}\t{'R1' if r1_match else 'R2'}\n")
            
            if total_reads >= max_reads:
                break

    print(f"Processed {total_reads} reads, matched {matched_reads} reads.")
    print(f"Output R1 FASTQ file: {out_r1_file}")
    print(f"Output R2 FASTQ file: {out_r2_file}")
    print(f"Output barcode file: {barcode_file}")

def main():
    parser = argparse.ArgumentParser(description='Process FASTQ files and extract reads with specific sequence.')
    parser.add_argument('-r1', '--read1', required=True, help='Input R1 FASTQ file (gzipped)')
    parser.add_argument('-r2', '--read2', required=True, help='Input R2 FASTQ file (gzipped)')
    parser.add_argument('-o', '--output', required=True, help='Output FASTQ file prefix')
    parser.add_argument('-m', '--max_reads', type=int, default=100000, 
                        help='Maximum number of total reads to process (default: %(default)s)')
    
    args = parser.parse_args()

    barcode_file = "Alpha.barcodes.txt"
    process_fastq_pair(args.read1, args.read2, args.output, barcode_file, args.max_reads)

if __name__ == "__main__":
    main()