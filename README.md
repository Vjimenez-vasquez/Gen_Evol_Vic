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

# codigo 6 #
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

# codigo 7 : Ensamblaje Nanopore (programas)
```r
# 7.1 : instalar los programas
# 7.1.1 : NanoPlot : calidad de secuencias Nanopore

conda install -c conda-forge -c bioconda nanoplot

or

pip install NanoPlot
pip install NanoPlot --upgrade

# 7.1.2 : Nanofilt : Filtrado por calidad de lecturas Nanopore

conda install -c bioconda nanofilt

or 

pip install nanofilt
pip install nanofilt --upgrade

# 7.1.3 : Flye: de-novo assembly

conda install -c bioconda flye

or 

git clone https://github.com/fenderglass/Flye
cd Flye
python setup.py install

# 7.1.4 : Minimap2 : polishing (parte 1)

conda install -c bioconda minimap2

or

git clone https://github.com/lh3/minimap2
cd minimap2 && make

# 7.1.5 : Racon : polishing (parte 2)

conda install -c bioconda racon

or 

git clone --recursive https://github.com/lbcb-sci/racon.git racon
cd racon
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make
cd build/bin/ 
export PATH=$PATH:$HOME/bin
cp racon $HOME/bin
chmod +x $HOME/bin/racon

# 7.1.6 : Requerimientos de MEDAKA (Pyabpoa, bcftools, samtools (v1.11), minimap2)

pip install pyabpoa
sudo apt install bcftools
conda install -c bioconda samtools==1.11

# 7.1.7 : MEDAKA, secuencias consenso (si MEDAKA no funciona correctamente, instala los programas requeridos)

conda install -c conda-forge –c bioconda medaka

or

pip install medaka
```

# codigo 8 : ensamblaje Nanopore (pipeline)
## material de apoyo > https://denbi-nanopore-training-course.readthedocs.io/en/stable/index.html ##
```r
# 8.1: descargar la informacion (códigos SRR17110067 y SRR17110070)
mkdir sra_files ;
prefetch --max-size 50G --option-file accessions.txt ;
mv */*.sra . ;
fasterq-dump --split-files *.sra 
gzip *.fastq ;
mkdir sra_files ;
mv *.sra sra_files/ ;

# 8.2: inspeccionar las longitudes de los reads ##
zcat SRR17110067.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20
zcat SRR17110070.fastq.gz | grep -n "length" | cut -f2 -d'=' | sort -r -n | uniq | head -n 20

# 8.3: NanoPlot
NanoPlot -t 2 -o SRR17110067_QC --fastq SRR17110067.fastq.gz
NanoPlot -t 2 -o SRR17110070_QC --fastq SRR17110070.fastq.gz

# 8.4: NanoFilt
gunzip -c SRR17110067.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110067.trim.fastq.gz ;
gunzip -c SRR17110070.fastq.gz | NanoFilt --logfile nanofilt.log -l 500 -q 10 | gzip > SRR17110070.trim.fastq.gz ;
ls -lh ;

# 8.5: Flye
flye -o SRR17110067.genoma --nano-raw SRR17110067.trim.fastq.gz --threads 4 ;
flye -o SRR17110070.genoma --nano-raw SRR17110070.trim.fastq.gz --threads 4 ;
ls -lh ;

# 8.6 : Minimap2 + Racon (Polishing)
minimap2 -x ava-ont -t 4 SRR17110067.genoma/assembly.fasta SRR17110067.trim.fastq.gz > overlaps1.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps1.paf SRR17110067.genoma/assembly.fasta > SRR17110067.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.genoma/assembly.fasta SRR17110070.trim.fastq.gz > overlaps2.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps2.paf SRR17110070.genoma/assembly.fasta > SRR17110070.racon1.fasta ;

minimap2 -x ava-ont -t 4 SRR17110067.racon1.fasta SRR17110067.trim.fastq.gz > overlaps3.paf ;
racon -t 4 SRR17110067.trim.fastq.gz overlaps3.paf SRR17110067.racon1.fasta > SRR17110067.racon2.fasta ;

minimap2 -x ava-ont -t 4 SRR17110070.racon1.fasta SRR17110070.trim.fastq.gz > overlaps4.paf ;
racon -t 4 SRR17110070.trim.fastq.gz overlaps4.paf SRR17110070.racon1.fasta > SRR17110070.racon2.fasta ;

# 8.7 : Medaka (consensus)
medaka_consensus -i SRR17110070.trim.fastq.gz -d SRR17110070.racon2.fasta -o medaka_SRR17110070 -t 4 ;
medaka_consensus -i SRR17110067.trim.fastq.gz -d SRR17110067.racon2.fasta -o medaka_SRR17110067 -t 4 ;

# 8.8 : QUAST
quast.py -o quast_results -m 0 consensus.fasta

# 8.9 : Bandage
```

# codigo 9 : BLAST
```r
# 9.1 : instalacion a traves de CONDA
conda install bioconda::blast

# 9.2 : http://www.mgc.ac.cn/VFs/
# 9.3 : Default webpage accessible to all users worldwide
# 9.4 : Download
# 9.5 : DNA sequences of full dataset
# 9.6 : Protein sequences of full dataset
# 9.6 :
gzip -d VFDB_setB_nt.fas.gz 
gzip -d VFDB_setB_pro.fas.gz

# 9.7 : run BLAST+
makeblastdb -in VFDB_setB_nt.fas -dbtype nucl ;
blastn -db VFDB_setB_nt.fas -query GCA_001183825.1.fasta -perc_identity 90 -outfmt 6 -num_threads 4 > blast.csv ;
head blast.csv ;
cat blast.csv ;

# 9.8 : headers
sed '1i query.acc.ver subject.acc.ver perc.identity alignment.length mismatches gap.opens q.start q.end s.start s.end evalue bit.score' blast.csv | tr " " "\t" > blast.2.csv

# 9.9 : revisar resultados
head blast.2.csv
cat blast.2.csv
```
