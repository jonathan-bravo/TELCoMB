#!/bin/bash

sample_filename=$1
amr_filename=$2
mges_filename=$3
raw_read_filename=$4
trim_filename=$5
hard_trim_filename=$6
chimeric_filename=$7
dedup_filename=$8
duplicate_filename=$9
non_host_read_filename=${10}
flow_cell=${11}
out=${12}

# Sample name
sample_name="${sample_filename##*/}"
sample_name="${sample_name%%_*}"
echo "$sample_name"

# AMR count
amr_count=$(grep -w "ARG_NUM_OF_READS" "$amr_filename" | cut -d, -f2 | tr -d '\r')
echo "$amr_count"

# MGE count
mges_count=$(grep "MGES_NUM_OF_READS" "$mges_filename" | cut -d, -f2 | tr -d '\r')
echo "$mges_count"

# Raw read count
raw_read_count=$(cat $raw_read_filename)
echo "$raw_read_count"

# Trim read count
trim_read_count=$(echo $(cat "$trim_filename" | wc -l) - 1 | bc)
echo "$trim_read_count"

# Hard trim read count
hard_trim_read_count=$(cat "$hard_trim_filename")
echo "$hard_trim_read_count"

# Chimeric read count
chimeric_read_count=$(cat "$chimeric_filename")
echo "$chimeric_read_count"

# Dedup read count
dedup_read_count=$(echo $(cat "$dedup_filename" | wc -l) - 1 | bc)
echo "$dedup_read_count"

# Duplicate read count
duplicate_read_count=$(echo $(cat "$duplicate_filename" | wc -l) | bc)
echo "$duplicate_read_count"

# Non-host read count
non_host_read_count=$(cat "$non_host_read_filename")
echo "$non_host_read_count"

# AMR on target raw reads
amr_ontarget_raw_read=$(echo "$amr_count $raw_read_count" | awk '{print $1 / $2}')
echo $amr_ontarget_rawread

# MGEs on target raw reads
mges_ontarget_raw_read=$(echo "$mges_count $raw_read_count" | awk '{print $1 / $2}')
echo $mges_ontarget_raw_read

# AMR on target non-host reads
amr_ontarget_non_host_read=$(echo "$amr_count $non_host_read_count" | awk '{print $1 / $2}')
echo $amr_ontarget_non_host_read

# MGEs on target non-host reads
mges_ontarget_non_host_read=$(echo "$mges_count $non_host_read_count" | awk '{print $1 / $2}')
echo $mges_ontarget_non_host_read

flow_cell_id="$flow_cell"
echo $flow_cell_id

# Concat counts
printf '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%f,%f,%f,%f\n' "$flow_cell_id" "$sample_name" "$amr_count" "$mges_count" "$raw_read_count" "$trim_read_count" "$hard_trim_read_count" "$chimeric_read_count" "$dedup_read_count" "$duplicate_read_count" "$non_host_read_count" "$amr_ontarget_raw_read" "$mges_ontarget_raw_read" "$amr_ontarget_non_host_read" "$mges_ontarget_non_host_read" > "$out"
