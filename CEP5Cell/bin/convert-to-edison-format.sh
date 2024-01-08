#!/bin/bash

# filename given from first argument
FILENAME=$1
# unit in ms or s
UNIT=$2
# number of column
COL="$(awk 'NR == 1{print NF}' $FILENAME)"
# number of rows
ROW="$(wc -l < $FILENAME)"
# header string which located in first line of original result file
HEADER="$(head -n 1 $FILENAME)"
# split HEADER string to array and save to HEADERS
IFS=' ' read -a HEADERS <<< $HEADER;

# print string number of field
echo "#NumField: $(($COL - 1))";
echo "#LabelX: Time ($UNIT), LabelY: Value";
# loop from second column to end
for i in $(seq 2 $COL); do
  # print field number, header label and number of rows
  echo "#Field$(($i - 1)): ${HEADERS[$(($i - 1))]},NumPoint:$(($ROW - 1))"
  # print content from original result file
  awk -v var="$i" 'NR > 1 {print $1 " " $var}' $FILENAME;
done
