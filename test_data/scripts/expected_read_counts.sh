#!/bin/bash

# Declare barcode array with corresponding read values
declare -A barcode_map

barcode_map["barcode01"]=$'Total reads\n10'
barcode_map["barcode02"]=$'Total reads\n2\n8'
barcode_map["barcode03"]=$'Total reads\n2\n8'
barcode_map["barcode04"]=$'Total reads\n12'
barcode_map["barcode05"]=$'Total reads\n2'
barcode_map["barcode06"]=$'Total reads\n2'
barcode_map["barcode07"]=$'Total reads\n2'
barcode_map["barcode08"]=$'Total reads\n2'
barcode_map["barcode09"]=$'Total reads\n2'
barcode_map["barcode10"]=$'Total reads\n1'
barcode_map["barcode11"]=$'Total reads\n1'
barcode_map["barcode12"]=$'Total reads\n1'

# Write each value to a separate file
for i in $(seq -w 1 12); do
    key="barcode$i"
    filename="test_data/${key}.tsv"
    echo -e "${barcode_map[$key]}" > "$filename"
done
