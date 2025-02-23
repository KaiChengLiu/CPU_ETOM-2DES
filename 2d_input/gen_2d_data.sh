#!/bin/bash

# 設定初始變數
t0=210
propagate_time=600

tau_step=10
tau_bound=600
input_file="key.key-tmpl"
T=(0)  # 使用 bash 陣列來存儲 T 的值

# 一次性讀取文件內容到變數中
input_content=$(cat "$input_file")

# 使用 awk 進行文本處理和替換
process_file() {
    local tau1=$1
    local tau2=$2
    local tau3=$3
    local t_end=$4
    local of_name=$5

    echo "Processing file: $of_name with tau1=$tau1, tau2=$tau2, tau3=$tau3, t_end=$t_end"

    # 用 awk 進行替換
    echo "$input_content" | awk -v tau1="$tau1" -v tau2="$tau2" -v tau3="$tau3" -v t_end="$t_end" '
    {
        gsub(/TAU1/, tau1);
        gsub(/TAU2/, tau2);
        gsub(/TAU3/, tau3);
        gsub(/T_END/, t_end);
        print;
    }' > "$of_name"
}

# 第一個迴圈：i 從 0 到 tau_bound
for t in "${T[@]}"; do
    for ((i=0; i<=tau_bound; i+=tau_step)); do
        tau1=$t0
        tau2=$(echo "$t0 + $i" | bc)
        tau3=$(echo "$t0 + $i + $t" | bc)
        t_end=$(echo "$t0 + $i + $t + $propagate_time" | bc)

        of_name="key_${i}_${t}.key"
        process_file "$tau1" "$tau2" "$tau3" "$t_end" "$of_name" &
    done
done

# 第二個迴圈：i 從 -tau_bound 到 0
for t in "${T[@]}"; do
    for ((i=-tau_bound; i<0; i+=tau_step)); do
        abs_i=$(echo "${i#-}")  # 取得 i 的絕對值
        tau1=$(echo "$t0 + $abs_i" | bc)
        tau2=$t0
        tau3=$(echo "$t0 + $abs_i + $t" | bc)
        t_end=$(echo "$t0 + $abs_i + $t + $propagate_time" | bc)

        of_name="key_${i}_${t}.key"
        process_file "$tau1" "$tau2" "$tau3" "$t_end" "$of_name" &
    done
done

# 等待所有背景進程完成
wait
