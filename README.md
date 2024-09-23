# Gen_Evol_Vic
una coleccion de codigos para la clase de Genomica Evolutiva

# codigo 1 : 
```r
#paso 1
https://github.com/ncbi/sra-tools/wiki/01.-Downloading-SRA-Toolkit
#paso 2
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz -O stk.tar.gz
#paso 3
chmod 777 stk.tar.gz
#paso 4
export PATH=$PATH:$PWD/sratoolkit.3.0.10-ubuntu64/bin
```

# codigo 2 
```r
#paso 1
prefetch --max-size 50G --option-file sra_accessions_1.txt ;
#paso 2 
prefetch -h
#paso 3
mv */*.sra . ;
#paso 4
rm -r ERR12389866/ ERR12543675/
#paso 5
fasterq-dump --split-files *.sra ;
#paso 6
gzip *fastq ;
#paso 7
fastqc *
```
