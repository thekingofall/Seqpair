import gzip
from Bio import SeqIO
import argparse
from tqdm import tqdm
import os

def count_reads(filename):
    with gzip.open(filename, 'rt') as f:
        return sum(1 for _ in f) // 4

def process_fastq_files(r1_file, r2_file, output_prefix, target_seq):
    # 计算总读段数
    total_reads = count_reads(r1_file)
    
    # 生成输出文件名
    output_r1 = f"{output_prefix}_R1.fastq.gz"
    output_r2 = f"{output_prefix}_R2.fastq.gz"
    
    # 打开输入和输出文件
    with gzip.open(r1_file, "rt") as handle1, \
         gzip.open(r2_file, "rt") as handle2, \
         gzip.open(output_r1, "wt") as out1, \
         gzip.open(output_r2, "wt") as out2:
        
        # 使用 SeqIO.parse 同时读取两个文件
        records1 = SeqIO.parse(handle1, "fastq")
        records2 = SeqIO.parse(handle2, "fastq")

        count = 0
        # 使用 tqdm 创建进度条
        for r1, r2 in tqdm(zip(records1, records2), total=total_reads, desc="Processing reads"):
            # 检查 R1 或 R2 是否包含目标序列
            if target_seq in r1.seq or target_seq in r2.seq:
                SeqIO.write(r1, out1, "fastq")
                SeqIO.write(r2, out2, "fastq")
                count += 1

    print(f"Total matching pairs: {count}")
    print(f"Output files: {output_r1} and {output_r2}")

def main():
    parser = argparse.ArgumentParser(description="Extract FASTQ reads containing a specific sequence.")
    parser.add_argument("--r1", required=True, help="Input R1 FASTQ.gz file")
    parser.add_argument("--r2", required=True, help="Input R2 FASTQ.gz file")
    parser.add_argument("--out", required=True, help="Output file prefix")
    parser.add_argument("--seq", default="ATATCCAGAACCCTGACCC", help="Target sequence to search for")
    args = parser.parse_args()

    process_fastq_files(args.r1, args.r2, args.out, args.seq)

if __name__ == "__main__":
    main()