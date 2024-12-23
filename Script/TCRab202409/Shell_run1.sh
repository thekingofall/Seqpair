#!/bin/bash

# 定义输入文件和输出文件
INPUT_FILE="TCRab_TCRab_R_combined_top10000.fq"
OUTPUT_FILE="sequence_analysis_results.csv"

# 函数：执行 grep 命令并返回结果数，处理简并碱基
count_occurrences() {
    local pattern="$1"
    # 将 M 替换为 [AC]，S 替换为 [GC]，K 替换为 [GT]
    pattern=$(echo "$pattern" | sed 's/M/[AC]/g' | sed 's/S/[GC]/g' | sed 's/K/[GT]/g')
    grep -ciP "$pattern" "$INPUT_FILE"
}

# 函数：计算多个序列组合的出现次数
count_multiple_occurrences() {
    local patterns=("$@")
    local grep_pattern=""
    for pattern in "${patterns[@]}"; do
        # 将 M 替换为 [AC]，S 替换为 [GC]，K 替换为 [GT]
        pattern=$(echo "$pattern" | sed 's/M/[AC]/g' | sed 's/S/[GC]/g' | sed 's/K/[GT]/g')
        grep_pattern+="(?=.*$pattern)"
    done
    grep -ciP "$grep_pattern" "$INPUT_FILE"
}

# 函数：生成反向互补序列
reverse_complement() {
    echo "$1" | tr 'ACGTMSKacgtmsk' 'TGCAKSMtgcaksm' | rev
}

# 定义序列
PT5A='GGGTCAGGGTTCTGGATAT'
P5TB='ACACSTTKTTCAGGTCCTC'
LINKER='AGCTGCTATGCACGACTG'

# 生成反向互补序列
PT5A_RC=$(reverse_complement "$PT5A")
P5TB_RC=$(reverse_complement "$P5TB")
LINKER_RC=$(reverse_complement "$LINKER")

# 计算单独序列出现次数
C_PT5A=$(count_occurrences "$PT5A")
C_P5TB=$(count_occurrences "$P5TB")
C_LINKER=$(count_occurrences "$LINKER")
C_PT5A_RC=$(count_occurrences "$PT5A_RC")
C_P5TB_RC=$(count_occurrences "$P5TB_RC")
C_LINKER_RC=$(count_occurrences "$LINKER_RC")

# 计算所有可能的组合
C_PT5A_LINKER=$(count_multiple_occurrences "$PT5A" "$LINKER")
C_PT5A_P5TB_RC=$(count_multiple_occurrences "$PT5A" "$P5TB_RC")
C_LINKER_P5TB_RC=$(count_multiple_occurrences "$LINKER" "$P5TB_RC")
C_PT5A_LINKER_P5TB_RC=$(count_multiple_occurrences "$PT5A" "$LINKER" "$P5TB_RC")

C_PT5A_RC_LINKER=$(count_multiple_occurrences "$PT5A_RC" "$LINKER")
C_PT5A_RC_P5TB=$(count_multiple_occurrences "$PT5A_RC" "$P5TB")
C_LINKER_P5TB=$(count_multiple_occurrences "$LINKER" "$P5TB")
C_PT5A_RC_LINKER_P5TB=$(count_multiple_occurrences "$PT5A_RC" "$LINKER" "$P5TB")

C_PT5A_RC_LINKER_RC=$(count_multiple_occurrences "$PT5A_RC" "$LINKER_RC")
C_LINKER_RC_P5TB=$(count_multiple_occurrences "$LINKER_RC" "$P5TB")
C_PT5A_RC_LINKER_RC_P5TB=$(count_multiple_occurrences "$PT5A_RC" "$LINKER_RC" "$P5TB")

C_P5TB_RC_LINKER_PT5A=$(count_multiple_occurrences "$P5TB_RC" "$LINKER" "$PT5A")

# 新增组合
C_PT5A_RC_LINKER_RC_P5TB_PT5A=$(count_multiple_occurrences "$PT5A_RC" "$LINKER_RC" "$P5TB" "$PT5A")
C_PT5A_RC_LINKER_RC_P5TB_P5TB_RC=$(count_multiple_occurrences "$PT5A_RC" "$LINKER_RC" "$P5TB" "$P5TB_RC")
C_PT5A_RC_LINKER_P5TB_P5TB_RC=$(count_multiple_occurrences "$PT5A_RC" "$LINKER" "$P5TB" "$P5TB_RC")
C_PT5A_LINKER_P5TB_P5TB_RC=$(count_multiple_occurrences "$PT5A" "$LINKER" "$P5TB" "$P5TB_RC")

# 创建 CSV 文件并写入表头
echo "序列组合,出现次数" > "$OUTPUT_FILE"

# 写入单独序列结果
echo "PT5A,$C_PT5A" >> "$OUTPUT_FILE"
echo "P5TB,$C_P5TB" >> "$OUTPUT_FILE"
echo "linker,$C_LINKER" >> "$OUTPUT_FILE"
echo "PT5A的反向互补,$C_PT5A_RC" >> "$OUTPUT_FILE"
echo "P5TB反向互补,$C_P5TB_RC" >> "$OUTPUT_FILE"
echo "linker反向互补,$C_LINKER_RC" >> "$OUTPUT_FILE"

# 添加分段符号
echo "---" >> "$OUTPUT_FILE"

# 写入组合序列结果
echo "PT5A + linker,$C_PT5A_LINKER" >> "$OUTPUT_FILE"
echo "PT5A + P5TB反向互补,$C_PT5A_P5TB_RC" >> "$OUTPUT_FILE"
echo "linker + P5TB反向互补,$C_LINKER_P5TB_RC" >> "$OUTPUT_FILE"
echo "PT5A + linker + P5TB反向互补,$C_PT5A_LINKER_P5TB_RC" >> "$OUTPUT_FILE"

echo "PT5A的反向互补 + linker,$C_PT5A_RC_LINKER" >> "$OUTPUT_FILE"
echo "PT5A的反向互补 + P5TB,$C_PT5A_RC_P5TB" >> "$OUTPUT_FILE"
echo "linker + P5TB,$C_LINKER_P5TB" >> "$OUTPUT_FILE"
echo "PT5A的反向互补 + linker + P5TB,$C_PT5A_RC_LINKER_P5TB" >> "$OUTPUT_FILE"

echo "PT5A的反向互补 + linker反向互补,$C_PT5A_RC_LINKER_RC" >> "$OUTPUT_FILE"
echo "linker反向互补 + P5TB,$C_LINKER_RC_P5TB" >> "$OUTPUT_FILE"
echo "**PT5A的反向互补 + linker反向互补 + P5TB,$C_PT5A_RC_LINKER_RC_P5TB" >> "$OUTPUT_FILE"

echo "**P5TB反向互补 + linker + PT5A,$C_P5TB_RC_LINKER_PT5A" >> "$OUTPUT_FILE"
echo "---" >> "$OUTPUT_FILE"

# 添加新的组合结果
echo "PT5A的反向互补 + linker反向互补 + P5TB + PT5A,$C_PT5A_RC_LINKER_RC_P5TB_PT5A" >> "$OUTPUT_FILE"
echo "PT5A的反向互补 + linker反向互补 + P5TB + P5TB的反向互补,$C_PT5A_RC_LINKER_RC_P5TB_P5TB_RC" >> "$OUTPUT_FILE"
echo "PT5A的反向互补 + linker + P5TB + P5TB的反向互补,$C_PT5A_RC_LINKER_P5TB_P5TB_RC" >> "$OUTPUT_FILE"
echo "PT5A + linker + P5TB + P5TB的反向互补,$C_PT5A_LINKER_P5TB_P5TB_RC" >> "$OUTPUT_FILE"

# 添加分段符号
echo "---" >> "$OUTPUT_FILE"

# 增加其他新的组合
C_PT5A_P5TB=$(count_multiple_occurrences "$PT5A" "$P5TB")
C_PT5A_LINKER_RC=$(count_multiple_occurrences "$PT5A" "$LINKER_RC")
C_P5TB_LINKER_RC=$(count_multiple_occurrences "$P5TB" "$LINKER_RC")
C_PT5A_P5TB_LINKER=$(count_multiple_occurrences "$PT5A" "$P5TB" "$LINKER")
C_PT5A_RC_P5TB_RC=$(count_multiple_occurrences "$PT5A_RC" "$P5TB_RC")

echo "PT5A + P5TB,$C_PT5A_P5TB" >> "$OUTPUT_FILE"
echo "PT5A + linker反向互补,$C_PT5A_LINKER_RC" >> "$OUTPUT_FILE"
echo "P5TB + linker反向互补,$C_P5TB_LINKER_RC" >> "$OUTPUT_FILE"
echo "PT5A + P5TB + linker,$C_PT5A_P5TB_LINKER" >> "$OUTPUT_FILE"
echo "PT5A的反向互补 + P5TB反向互补,$C_PT5A_RC_P5TB_RC" >> "$OUTPUT_FILE"

echo "分析完成。结果已保存到 $OUTPUT_FILE"