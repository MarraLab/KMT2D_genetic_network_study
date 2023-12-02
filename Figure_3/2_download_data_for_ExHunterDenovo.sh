#############################################################
# Figure 3
# Extract SRA fastq files for Expansion Hunter Denovo analysis
#############################################################
# # Sample SRA ID lists
all_wgs_samples=(SRR8639191 SRR8652098 SRR8670708 SRR8652101 SRR8652121 SRR8652136 SRR8670772 SRR11680468 SRR11680467 SRR8670762 SRR8670679 SRR8670673 SRR8670700 SRR8639227 SRR8639145 SRR8652076)

# Prefetch 
for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists."
        /projects/marralab/ytakemon_prj/Programs/sratoolkit.2.11.3-centos_linux64/bin/prefetch.2.11.2 --max-size 200G -p $SRA_ID > prefetch.log 2>&1 &
    fi
done

# Convert .sra file to .fastq the rest of the dataset 
for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists."
        fasterq-dump $SRA_ID -pv --threads 10 -t temp_sra_dump > fasterq-dump.log 2>&1 &
    fi
done

# Compress fastq
for file in *.fastq
do 
    gzip ${file} &
done

# Align fastq
cd 
hg38=alignment_references/Homo_sapiens/hg38_no_alt/genome/bwa_64/hg38_no_alt.fa
bwa=linux-x86_64-centos7/bwa-0.7.17/bwa
samtools=linux-x86_64-centos7/samtools-1.9/bin/samtools

# Testing bash statemets 
for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists."
    fi
done

# align, convert to bam, sort, and index. 
for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists. Creating SAM file"
        ${bwa} mem -t9 ${hg38} ${SRA_ID}_1.fastq.gz ${SRA_ID}_2.fastq.gz > ${SRA_ID}.sam &
    fi    
done

for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists. Creating BAM file"
        samtools view -@9 -b ${SRA_ID}.sam -o ${SRA_ID}.bam &
    fi    
done

for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists. Creating sorted BAM file"
        samtools sort -@9 ${SRA_ID}.bam > ${SRA_ID}.sorted.bam &
    fi    
done


for SRA_ID in ${all_wgs_samples[@]}
do
    FILE=(${SRA_ID}.sorted.bam.bai)
    if [ -f "$FILE" ]; then
        echo "$FILE exists."
    else 
        echo "$FILE doesn't exists. Creating index file"
        samtools index -@9 ${SRA_ID}.sorted.bam &
    fi    
done

# clean up
for SRA_ID in ${all_wgs_samples[@]}
do
    rm ${SRA_ID}.fastq.gz &
    rm ${SRA_ID}.sam &
    rm ${SRA_ID}.bam &