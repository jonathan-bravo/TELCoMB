#!/bin/bash

sample_filename=$1
amr_filename=$2
mges_filename=$3
raw_read_filename=$4
non_host_read_filename=$5
out=$6

sample_name="${sample_filename##*/}"
sample_name="${sample_name%%_*}"
echo "$sample_name"

amr_count=$(grep -w "ARG_NUM_OF_READS" "$amr_filename" | cut -d, -f2 | tr -d '\r')
echo "$amr_count"

mges_count=$(grep "MGES_NUM_OF_READS" "$mges_filename" | cut -d, -f2 | tr -d '\r')
echo "$mges_count"

raw_read_count=$(cat $raw_read_filename)
echo "$raw_read_count"

non_host_read_count=$(cat $non_host_read_filename)
echo "$non_host_read_count"

amr_ontarget_raw_read=$(echo "$amr_count $raw_read_count" | awk '{print $1 / $2}')
echo $amr_ontarget_rawread

mges_ontarget_raw_read=$(echo "$mges_count $raw_read_count" | awk '{print $1 / $2}')
echo $mges_ontarget_raw_read

amr_ontarget_non_host_read=$(echo "$amr_count $non_host_read_count" | awk '{print $1 / $2}')
echo $amr_ontarget_non_host_read

mges_ontarget_non_host_read=$(echo "$mges_count $non_host_read_count" | awk '{print $1 / $2}')
echo $mges_ontarget_non_host_read

printf '%s,%s,%s,%s,%s,%f,%f,%f,%f\n' "$sample_name" "$amr_count" "$mges_count" "$raw_read_count" "$non_host_read_count" "$amr_ontarget_raw_read" "$mges_ontarget_raw_read" "$amr_ontarget_non_host_read" "$mges_ontarget_non_host_read" > "$out"
