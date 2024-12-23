# mkdir -p Top10000 && for file in *_1.fq.gz; do base=$(basename $file _1.fq.gz); paste <(zcat ${base}_1.fq.gz | head -n 40000) <(zcat ${base}_2.fq.gz | head -n 40000) > Top10000/${base}_combined_top10000.fq; done
import os
import regex
import argparse
import pandas as pd

def count_lines(filepath):
    with open(filepath, 'rt') as f:
        return sum(1 for _ in f)

def process_file(fastq_path, search_sequence, tolerance):
    # 定义搜索的序列，将S和M的模糊匹配转换为正则表达式
    search_sequence = search_sequence.replace('S', '[CG]').replace('M', '[AC]')

    # 初始化计数器
    sequence_counts = {
        'forward': 0,
        'complement': 0,
        'reverse_complement': 0,
        'reverse': 0,
        'non_matching': 0
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

    # 获取文件总行数用于总读取数计算
    total_lines = count_lines(fastq_path)
    total_reads = total_lines // 4  # 每个读取由4行组成

    # 读取和处理 FASTQ 文件
    with open(fastq_path, 'rt') as fastq:
        for line_number, line in enumerate(fastq):
            if line_number % 4 == 1:  # 序列行
                sequence = line.strip()
                matched = False
                for pattern_name, pattern in regex_patterns.items():
                    if pattern.search(sequence):
                        sequence_counts[pattern_name] += 1
                        matched = True
                        break
                if not matched:
                    sequence_counts['non_matching'] += 1

    # 生成输出内容
    result = {
        'search_sequence': search_sequence,
        'filename': os.path.basename(fastq_path),
        'total_reads': total_reads
    }
    for seq_type, count in sequence_counts.items():
        percentage = (count / total_reads) * 100
        result[seq_type] = f"{count} ({percentage:.2f}%)"

    return result

def main(input_folder, output_file, search_sequence, tolerance):
    # 获取目录中所有 .fq 文件
    fastq_files = [os.path.join(input_folder, f) for f in os.listdir(input_folder) if f.endswith('.fq')]

    # 初始化结果列表
    results = []

    # 处理每个文件
    for fastq_file in fastq_files:
        print(f"Processing {fastq_file}...")
        result = process_file(fastq_file, search_sequence, tolerance)
        results.append(result)

    # 转换为 Pandas DataFrame
    df = pd.DataFrame(results)

    # 打印到控制台
    print(df)

    # 写入结果到输出文件
    df.to_csv(output_file, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Count sequences in FASTQ files within a folder.')
    parser.add_argument('-i', '--input_folder', type=str, required=True, help='Path to the folder containing FASTQ files.')
    parser.add_argument('-o', '--output_file', type=str, required=True, help='Output CSV file to write the results.')
    parser.add_argument('-s', '--search_sequence', type=str, default='GAGGACCTGA[AC][AC][CG]A[CG]TGT', help='Search sequence with S and M ambiguity codes.')
    parser.add_argument('-t', '--tolerance', type=int, default=2, help='Number of tolerated mismatches.')

    args = parser.parse_args()

    main(args.input_folder, args.output_file, args.search_sequence, args.tolerance)