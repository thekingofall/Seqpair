import gzip
import os
import shutil
import argparse
from Bio import SeqIO
from tqdm import tqdm
from multiprocessing import Pool, cpu_count
import itertools

def count_reads(filename):
    with gzip.open(filename, 'rt') as f:
        return sum(1 for _ in f) // 4

def process_chunk(args):
    chunk, r1_file, r2_file, temp_dir, target_seq, chunk_size, chunk_number = args
    output_r1 = os.path.join(temp_dir, f"chunk{chunk_number}_R1.fastq.gz")
    output_r2 = os.path.join(temp_dir, f"chunk{chunk_number}_R2.fastq.gz")
    
    count = 0
    with gzip.open(r1_file, "rt") as handle1, \
         gzip.open(r2_file, "rt") as handle2, \
         gzip.open(output_r1, "wt") as out1, \
         gzip.open(output_r2, "wt") as out2:
        
        records1 = itertools.islice(SeqIO.parse(handle1, "fastq"), chunk * chunk_size, (chunk + 1) * chunk_size)
        records2 = itertools.islice(SeqIO.parse(handle2, "fastq"), chunk * chunk_size, (chunk + 1) * chunk_size)

        for r1, r2 in zip(records1, records2):
            if target_seq in r1.seq or target_seq in r2.seq:
                SeqIO.write(r1, out1, "fastq")
                SeqIO.write(r2, out2, "fastq")
                count += 1

    return count, output_r1, output_r2

def process_fastq_files(r1_file, r2_file, output_prefix, target_seq, num_threads, temp_dir):
    total_reads = count_reads(r1_file)
    chunk_size = total_reads // num_threads + 1
    
    # 创建临时目录
    os.makedirs(temp_dir, exist_ok=True)
    
    pool = Pool(num_threads)
    chunks = range(num_threads)
    args = [(chunk, r1_file, r2_file, temp_dir, target_seq, chunk_size, chunk) for chunk in chunks]
    
    total_count = 0
    output_files_r1 = []
    output_files_r2 = []
    
    with tqdm(total=total_reads, desc="Processing reads") as pbar:
        for result in pool.imap_unordered(process_chunk, args):
            count, out_r1, out_r2 = result
            total_count += count
            output_files_r1.append(out_r1)
            output_files_r2.append(out_r2)
            pbar.update(chunk_size)
    
    pool.close()
    pool.join()

    # Merge output files
    final_output_r1 = f"{output_prefix}_R1.fastq.gz"
    final_output_r2 = f"{output_prefix}_R2.fastq.gz"
    
    print("Merging output files...")
    with gzip.open(final_output_r1, 'wb') as outfile:
        for fname in output_files_r1:
            with gzip.open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)
    
    with gzip.open(final_output_r2, 'wb') as outfile:
        for fname in output_files_r2:
            with gzip.open(fname, 'rb') as infile:
                shutil.copyfileobj(infile, outfile)

    print(f"Total matching pairs: {total_count}")
    print(f"Output files: {final_output_r1} and {final_output_r2}")

    # 可选：删除临时文件夹
    shutil.rmtree(temp_dir)

def main():
    parser = argparse.ArgumentParser(description="Extract FASTQ reads containing a specific sequence.")
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ.gz file")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ.gz file")
    parser.add_argument("--out", required=True, help="Output file prefix")
    parser.add_argument("--seq", default="ATATCCAGAACCCTGACCC", help="Target sequence to search for")
    parser.add_argument("--threads", type=int, default=cpu_count(), help="Number of threads to use")
    parser.add_argument("--temp", default="temp_fastq", help="Temporary directory for intermediate files")
    args = parser.parse_args()

    process_fastq_files(args.r1, args.r2, args.out, args.seq, args.threads, args.temp)

if __name__ == "__main__":
    main()