#!/bin/bash

genome=$1
for file in $(cat $genome)
do 
datasets download genome accession ${file} --filename ${file}.zip ;
unzip ${file}.zip -d ${file}/ ;
rm ${file}.zip ;
cp ${file}/*/*/*/*.fna .
rm -r ${file}/ ;
done

for file in *.fna
do
mv "$file" "$(echo "$file" | sed -e 's/.1_.*/.1.fasta/g')"
mv "$file" "$(echo "$file" | sed -e 's/.2_.*/.2.fasta/g')"
mv "$file" "$(echo "$file" | sed -e 's/.3_.*/.3.fasta/g')"
done 
ls ;
