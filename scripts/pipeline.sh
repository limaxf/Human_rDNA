module purge
module load sra-tools/3.1.0

#reconfigure the sratool target directory
vdb-config --prefetch-to-cwd

module load bowtie2/2.4.4
module load lofreq/2.1.5

#this module swap is important to run before loading samtools(suspect it has something to do with the multiple versions available
module swap htslib/intel/1.12 htslib/intel/1.14

module load samtools/intel/1.14
module load bedtools/intel/2.29.2

###########################################################################

# Loop over the sample identifiers
filename="$1"

echo "$filename"
#read -r -a sample_list < "$filename"
sample_list=()
while IFS= read -r line; do
  sample_list+=("$line")
done < "$filename"
echo "$sample_list"

for sample in "${sample_list[@]}"
do
    start1=$(date +%s.%N)
    echo "Processing sample: $sample"
#    if [ -d $sample ]; then
#	echo "Folder 'sample' already exists."
#    else
#        echo "Folder 'sample' does not exist."
#    fi
    #fetchNdump
    prefetch ${sample}
    fasterq-dump ${sample}
    rm -r ${sample}

    #align reads to rDNA
    bowtie2 -5 1 -N 1 -p 8 \
    -x /scratch/cgsb/hochwagen/Human_rDNA_project/rDNA_prototype_prerRNA_only/rDNA_prototype_prerRNA_only \
    -1 ${sample}_1.fastq \
    -2 ${sample}_2.fastq \
    -S ${sample}_output.sam
    rm ${sample}_1.fastq
    rm ${sample}_2.fastq

    #sort
    samtools view -Sbh ${sample}_output.sam > ${sample}_rDNA.bam
    samtools sort -@ 8 -o ${sample}_sort.bam -O 'bam' ${sample}_rDNA.bam
    rm ${sample}_output.sam
    rm ${sample}_rDNA.bam


    ## convert to cram format
    samtools view -C -T /scratch/cgsb/hochwagen/Human_rDNA_project/rDNA_prototype_prerRNA_only.fa \
    -o ${sample}_rDNA.cram ${sample}_sort.bam

    #lofreq
    lofreq indelqual --dindel \
    -f /scratch/cgsb/hochwagen/Human_rDNA_project/rDNA_prototype_prerRNA_only.fa \
    -o ${sample}_rDNA.bam \
    ${sample}_rDNA.cram
    
  
    #calculate coverage
    bedtools genomecov -d \
    -ibam ${sample}_rDNA.bam > ${sample}_rDNA_coverage.txt

    #call variants
    lofreq call --call-indels -f /scratch/cgsb/hochwagen/Human_rDNA_project/rDNA_prototype_prerRNA_only.fa \
    -o ${sample}_rDNA.vcf \
    ${sample}_rDNA.bam 
    

    rm ${sample}_sort.bam
    rm ${sample}_rDNA.bam
    rm ${sample}_rDNA.cram 
####move files to it corresponding folder.
    
    mkdir ${sample}
    mv ${sample}_rDNA.vcf ${sample}/
#   mv ${sample}_rDNA.bam ${sample}/
    #mv ${sample}_rDNA.cram ${sample}/
    mv ${sample}_rDNA_coverage.txt ${sample}/
    echo "Processing of sample $sample complete"
    dur1=$(echo "$(date +%s.%N) - $start1" | bc)
    printf "Execution time: %.6f seconds" $dur1
done
