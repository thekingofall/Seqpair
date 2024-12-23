import gzip
import os
import regex
import sys
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from tqdm import tqdm
import glob

# Function to create a regex pattern allowing for ambiguous bases and mismatches
def generate_pattern(sequence, max_mismatches):
    sequence = (sequence.replace('S', '[CG]')
                        .replace('M', '[AC]'))
    # Allow up to max_mismatches
    return regex.compile(f"({sequence}){{e<={max_mismatches}}}")

# Function to open the fastq.gz files
def open_fastq(path):
    return gzip.open(path, 'rt')

# Function to write to fastq.gz files
def write_fastq(fastq_path, records):
    with gzip.open(fastq_path, 'wt') as f:
        for title, sequence, quality in records:
            f.write(f"@{title}\n{sequence}\n+\n{quality}\n")

# Function to process reads and output modified files
# Function to process reads and output modified files
def process_reads(input_basename, search_sequence, max_mismatches, output_folder):
    # Generate the pattern for reverse complement search
    rev_complement_pattern = generate_pattern(str(Seq(search_sequence).reverse_complement()), max_mismatches)

    # Lists to hold modified records
    modified_records_r1 = []
    modified_records_r2 = []

    # Count processed pairs
    pair_count = 0

    # Open both R1 and R2 files
    with open_fastq(f"{input_basename}_R1.fastq.gz") as fastq_r1, open_fastq(f"{input_basename}_R2.fastq.gz") as fastq_r2:
        for (title_r1, sequence_r1, quality_r1), (title_r2, sequence_r2, quality_r2) in tqdm(zip(FastqGeneralIterator(fastq_r1), FastqGeneralIterator(fastq_r2))):
            # Search for the reverse complement pattern in both reads
            match_r1 = rev_complement_pattern.search(sequence_r1)
            match_r2 = rev_complement_pattern.search(sequence_r2)

            # If either read contains the search sequence, trim and add to their respective modified records
            if match_r1 or match_r2:
                if match_r1:
                    end_r1 = match_r1.end()  # Get the end position of the match in R1
                    sequence_r1 = sequence_r1[end_r1:]  # Trim the sequence
                    quality_r1 = quality_r1[end_r1:]  # Trim the quality scores
                if match_r2:
                    end_r2 = match_r2.end()  # Get the end position of the match in R2
                    sequence_r2 = sequence_r2[end_r2:]  # Trim the sequence
                    quality_r2 = quality_r2[end_r2:]  # Trim the quality scores

                # Add the modified records to the lists
                modified_records_r1.append((title_r1, sequence_r1, quality_r1))
                modified_records_r2.append((title_r2, sequence_r2, quality_r2))

            # Increment the pair counter
            # pair_count += 1
            # #Stop after processing 10000 read pairs
            # if pair_count >= 300000:
            #     break 

    # Write the modified records to new fastq.gz files
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    trimmed_r1_path = os.path.join(output_folder, "trimmed_1.fastq.gz")
    trimmed_r2_path = os.path.join(output_folder, "trimmed_2.fastq.gz")
    write_fastq(trimmed_r1_path, modified_records_r1)
    write_fastq(trimmed_r2_path, modified_records_r2)

    # Call the TCR_get function
    TCR_get(output_folder, output_folder)

# Placeholder for the TCR_get function

def TCR_get(readFolder,outdirname):
        # TCRdir=os.path.join(outdirname,"TCR")
        TCRruntxt=os.path.join(outdirname,"TCRrun.txt")
        TCRrunt=open(TCRruntxt,"w")
        for file in glob.glob(readFolder+"/*_1.f*"):
                samplename=file.split("/")[-1].split("_1.f")[0]
                print(samplename)

                read1=file.split("/")[-1]
                # print(samplename)
                print(file.split("/")[-1].split("_")[0])
                read2=file.split("/")[-1].split("_1.f")[0]+"_2.f"+file.split("/")[-1].split("_1.f")[1]
                read2=read2.replace("_1_","_2_")
                # print(read2.replace("_1_","_2_"))
                print(read2)
                if os.path.exists(readFolder+"/"+read2) and read2.endswith(".gz") and read1.endswith(".gz"):
                        read2=read2
                        read1=read1
                        outdir=os.path.join(outdirname,samplename+"_TCR")
                        outpre=samplename+"_TCRtype"
                        
                        os.system("mkdir -p "+outdir)
                        TCRruncommnd="java -jar /home/maolp/mao/Project/20211110_xcc_TCR/2.cleandata/mixcr-3.0.13/mixcr.jar  analyze shotgun -s hsa --starting-material RNA  -t 20  "+readFolder+"/"+read1+" "+readFolder+"/"+read2 +" "+outdir+"/"+outpre
                        os.system(TCRruncommnd)
                        print(TCRruncommnd)
                        TCRrunt.write(TCRruncommnd+"\n")

# Main function to run the script
def main():
    if len(sys.argv) != 5:
        print("Usage: python script.py <input_fastq_basename> <sequence_to_search> <max_mismatches> <output_folder>")
        sys.exit(1)

    input_basename = sys.argv[1]
    search_sequence = sys.argv[2].upper()
    max_mismatches = int(sys.argv[3])
    output_folder = sys.argv[4]

    # Call the processing function
    process_reads(input_basename, search_sequence, max_mismatches, output_folder)

# Run the main function
if __name__ == "__main__":
    main()