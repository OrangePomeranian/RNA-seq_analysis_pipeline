#!/bin/bash

### Step 1: Download sequencing data
echo "Fetching sequencing run information..."
esearch -db sra -query PRJNA313294 | efetch -format runinfo -mode xml | xtract -pattern SraRunInfo -element Run > info.txt

echo "Processing run information..."
awk -F'\t' '{for(i=1;i<=NF;i++) print $i}' info.txt > runinfo.txt

echo "Downloading FASTQ files..."
while read p; do
    AR=${p}
    echo "Downloading: ${AR}"
    fastq-dump ${AR} --split-files
done < runinfo.txt

### Step 2: Quality control with FastQC
echo "Running FastQC on downloaded FASTQ files..."
fastqc *.fastq

echo "Organizing FastQC reports..."
mkdir -p qcfast
mv *.html qcfast/
mv *.zip qcfast/

### Step 3: Trimming sequences with Trimmomatic
echo "Trimming sequences for quality control..."

# Define arrays for single-end and paired-end reads
single=("SRR3194428_1" "SRR3194429_1" "SRR3194430_1" "SRR3194431_1")
pair1=("SRR3191542_1" "SRR3191543_1" "SRR3191544_1" "SRR3191545_1")
pair2=("SRR3191542_2" "SRR3191543_2" "SRR3191544_2" "SRR3191545_2")

echo "Trimming single-end reads..."
for file in "${single[@]}"; do
    trimmomatic SE "$file.fastq" "$file.fq" SLIDINGWINDOW:4:30 TRAILING:30
done

echo "Trimming paired-end reads..."
for index in "${!pair1[@]}"; do
    r1=${pair1[$index]}
    r2=${pair2[$index]}
    trimmomatic PE "$r1.fastq" "$r2.fastq" "$r1.fq" "${r1}_unpaired.fq" "$r2.fq" "${r2}_unpaired.fq" SLIDINGWINDOW:4:30 TRAILING:30
done

echo "Organizing pre-trimming FASTQ files..."
mkdir -p pretrimmomatic
mv *.fastq pretrimmomatic/

### Step 4: Indexing reference genome for HISAT2
echo "Building HISAT2 index..."
hisat2-build -f -p 6 hg19.fa hg19

### Step 5: Aligning reads using HISAT2
echo "Aligning single-end reads..."
for file in "${single[@]}"; do
    echo "Aligning $file"
    hisat2 -p 6 -x hg19 "$file.fq" -S "$file.single.sam"
    echo "$file alignment done"
done

echo "Aligning paired-end reads..."
for index in "${!pair1[@]}"; do
    r1=${pair1[$index]}
    r2=${pair2[$index]}
    echo "Aligning $r1 and $r2"
    hisat2 -q -p 6 -x hg19 -1 "$r1.fq" -2 "$r2.fq" -S "$r1.sam"
    echo "$r1 alignment done"
done

### Step 6: Convert SAM to BAM and sort BAM files
listasambam=("SRR3194428_1.single" "SRR3194429_1.single" "SRR3194430_1.single" "SRR3194431_1.single" "SRR31
91542_1" "SRR3191543_1" "SRR3191544_1" "SRR3191545_1")

echo "Converting SAM to BAM..."
for file in "${listasambam[@]}"; do
    echo "Processing $file"
    samtools view -Sb -@ 6 "$file.sam" > "$file.bam"
done

echo "Sorting BAM files..."
for file in "${listasambam[@]}"; do
    samtools sort "$file.bam" -o "${file}sorted.bam"
done

### Step 7: Counting features using featureCounts
echo "Generating count matrix..."
featureCounts -a hg19.gtf -o counts.txt *sorted.bam
featureCounts -a hg19.gtf -o countsSe.txt *singlesorted.bam 

echo "Pipeline completed successfully!"
