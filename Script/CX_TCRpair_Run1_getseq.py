import gzip
import regex
import argparse
from tqdm import tqdm

def count_lines(filepath):
    with gzip.open(filepath, 'rt') as f:
        return sum(1 for _ in f)

def main(fastq_gz_path, output_path, search_sequence, tolerance):
    # 定义搜索的序列，将S和M的模糊匹配转换为正则表达式
    search_sequence = search_sequence.replace('S', '[CG]').replace('M', '[AC]')

    # 初始化计数器
    sequence_counts = {
        'forward': 0,
        'complement': 0,
        'reverse_complement': 0,
        'reverse': 0
    }

    # 计算互补和反向互补模式
    complement = str.maketrans('ACGT', 'TGCA')
    forward_pattern = search_sequence
    reverse_pattern = forward_pattern[::-1]
    complement_pattern = forward_pattern.translate(complement)
    reverse_complement_pattern = reverse_pattern.translate(complement)

    # 生成正则表达式
    fuzziness = f'{{e<={tolerance}}}'
    regex_patterns = {
        'forward': regex.compile(f'({forward_pattern}){fuzziness}'),
        'reverse': regex.compile(f'({reverse_pattern}){fuzziness}'),
        'complement': regex.compile(f'({complement_pattern}){fuzziness}'),
        'reverse_complement': regex.compile(f'({reverse_complement_pattern}){fuzziness}')
    }

    # 获取文件总行数用于进度条
    total_lines = count_lines(fastq_gz_path)

    # 读取和处理 FASTQ 文件
    with gzip.open(fastq_gz_path, 'rt') as fastq:
        for line_number, line in tqdm(enumerate(fastq), total=total_lines, desc='Processing', unit='line'):
            if line_number % 4 == 1:  # 序列行
                sequence = line.strip()
                for pattern_name, pattern in regex_patterns.items():
                    if pattern.search(sequence):
                        sequence_counts[pattern_name] += 1

    # 写入结果到输出文件
    with open(output_path, 'w') as output_file:
        output_file.write("Sequence counts in FASTQ file:\n")
        for seq_type, count in sequence_counts.items():
            output_file.write(f"{seq_type.capitalize()} sequence count: {count}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count sequences in FASTQ file.')
    parser.add_argument('fastq_gz_path', type=str, help='Path to the FASTQ.GZ file.')
    parser.add_argument('output_path', type=str, help='Output file to write the results.')
    parser.add_argument('-s', '--search_sequence', type=str, default='GAGGACCTGA[AC][AC][CG]A[CG]TGT', help='Search sequence with S and M ambiguity codes.')
    parser.add_argument('-t', '--tolerance', type=int, default=2, help='Number of tolerated mismatches.')

    args = parser.parse_args()

    main(args.fastq_gz_path, args.output_path, args.search_sequence, args.tolerance)