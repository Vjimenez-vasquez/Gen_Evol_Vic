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
# codigo 3 
```r
#1# indexar el genoma de referencia#
bwa index reference.fasta ;

#2# preparar las instrucciones generales#
for r1 in *fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz

#3# instrucciones para generar el archivo .bam#
bwa mem -t 4 reference.fasta $r1 $r2 > ${prefix}_uno.sam ;
samtools view -@ 4 -bS -T reference.fasta ${prefix}_uno.sam > ${prefix}_unoa.bam ;
samtools sort -@ 4 -n ${prefix}_unoa.bam -o ${prefix}_dosa.bam ;
samtools fixmate -@ 4 -m ${prefix}_dosa.bam ${prefix}_tresa.bam ;
samtools sort -@ 4 ${prefix}_tresa.bam -o ${prefix}_cuatroa.bam ;
samtools markdup -@ 4 ${prefix}_cuatroa.bam ${prefix}.bam ;
samtools index -@ 4 ${prefix}.bam ;
rm ${prefix}_uno.sam ${prefix}_unoa.bam ${prefix}_dosa.bam ${prefix}_tresa.bam ${prefix}_cuatroa.bam ;
done ;
ls ;

#4# extraer genoams consenso #
for r1 in *bam
do
prefix=$(basename $r1 .bam)
#2#estimate Ns#
samtools mpileup -aa -A -d 0 -Q 0 $r1 | ivar consensus -p ${prefix}.fasta -q 25 -t 0.6 -m 10 ;
done ; 
ls ;
```

# codigo 4
```r
# instalacion de prokka #
conda create -n prokka_env ;
conda activate prokka_env ;
conda install -c conda-forge -c bioconda prokka ;

# analisis en prokka #
mkdir annotation/ ;
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;

# instalacion de artemis #
conda create -n art
conda activate art
conda install bioconda::artemis
```

# codigo 5 #
```r
# si desea correr el comando desde un archico ".sh" debe añardir la "shebang" al inicio #
# 5.1: crear un archivo con extension *.sh, ejemplo "comando_1.sh"#
# 5.2: pegar el siguiente contenido o a través de cat o a través de "nano" #

#!/usr/bin
for r1 in *fa
do
prefix=$(basename $r1 .fa)
prokka --cpus 4 $r1 -o ${prefix} --prefix ${prefix} --kingdom Viruses ; 
mv ${prefix}/*.gff annotation/${prefix}.gff
done ;

# 5.3: dar permiso al archivo generado
chmod 777 comando_1.sh

# 5.4: correr el programa
./comando_1.sh
```

## codigo 6 ##
```r
# 6.1: para obervar los headers de cada contig en todos los archivos *.fa
grep ">" *.fa

# 6.2: para observar todo el contenido de todos los archivos *.fa
cat *.fa

# 6.3: para obervar las 10 primeras lineas de cada archivo *.fa
cat OQ603638.fa | head -n 10
cat OQ603651.fa | head -n 10
cat SRR30716253.fa | head -n 10
cat SRR30716253.fa | head -n 10

# 6.4: para observar los resultados de la anotación en ARTEMIS, debe contar con el archivo *fa original y el  archivo *.gff. Cargue primero el genoma y luego el archivo de anotación
conda activate art
art
file >> open file manager >> cargar el genoma en extension *fa
file >> read and entry >> cargar el archivo *gff
explorar

# 6.5: identifique las regiones inferidas por el programa (CDS), identifique si se identificó la identidad de esas regiones o si algunas aparecen como "hypothetical"
```
