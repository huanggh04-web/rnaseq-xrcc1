#!/usr/bin/env bash
# 统计每个样本中「完美匹配人类 XRCC1」的 reads 数（= 转基因表达量）
cd ~/projects/rnaseq_xrcc1
mkdir -p results/human_xrcc1
NREADS=16000000        # 每个样本取前 400 万条 R1 reads

count_one() {
  local f="$1" s
  s=$(basename "$f" | sed 's/\..*//')
  zcat "$f" | head -n "$NREADS" \
    | bowtie2 -x refs/human/hXRCC1 -U - --very-sensitive -p 4 \
        --no-hd --no-sq 2>/dev/null \
    | awk '$1 !~ /^@/ {
        n++
        nm = -1
        for (i = 12; i <= NF; i++) if ($i ~ /^NM:i:/) { nm = substr($i,6)+0; break }
        if (nm < 0) next
        a++
        if      (nm == 0)  b0++
        else if (nm <= 3)  b3++
        else               bhi++
      }
      END { printf "%d\t%d\t%d\t%d\t%d\n", n+0, a+0, b0+0, b3+0, bhi+0 }' \
    > "results/human_xrcc1/$s.tsv"
  echo "  done: $s"
}
export -f count_one
export NREADS

ls ~/rawdata/*R1.raw.fastq.gz | sort \
  | xargs -P 6 -I{} bash -c 'count_one "$@"' _ {}
