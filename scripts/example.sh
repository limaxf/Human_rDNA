########################################################################
Reminder: watch the available space, files could be huge
########################################################################

####1. filter reads that are in the rDNA region
samtools view -h wgs_1011086.cram chr21:8433222-8446571 chr21:8205988-8219301 chr21:8389035-8402342 chr22_KI270733v1_random:122273-135644 chrUn_GL000220v1:105424-118779 > RNA45SN1t5.cram
####2. reconstruct reads in paried fastq format
samtools collate -u -O RNA45SN1t5.cram | \
samtools fastq -1 paired1.fq -2 paired2.fq -0 /dev/null -s /dev/null -n
####3. alignment to our reference
/home/jupyter/packages/bowtie2-2.5.4-linux-x86_64/bowtie2 -5 1 -N 1 -p 8 \
    -x /home/jupyter/workspaces/humanrrnaproject/rDNA_prototype_prerRNA/rDNA_prototype_prerRNA_only \
    -1 paired1.fq \
    -2 paired2.fq \
    -S example_output.sam
####4. convert to bam
samtools view -Sbh example_output.sam > example_output_rDNA.bam
####5.miniconda activation for lofreq
source ~/miniconda3/bin/activate 
conda activate lofreq
####6. mark indel
lofreq indelqual --dindel \
    -f /home/jupyter/workspaces/humanrrnaproject/rDNA_prototype_prerRNA_only.fa \
    -o example_indel_rDNA.bam \
    example_output_rDNA.bam
####7. sort (sort at the last step, lofreq needs it sorted)
samtools sort -@ 8 -o example_indel_sort.bam -O 'bam'  example_indel_rDNA.bam
####8. call variant
lofreq call --call-indels -f /home/jupyter/workspaces/humanrrnaproject/rDNA_prototype_prerRNA_only.fa \
    -o example_output_rDNA.vcf \
    example_indel_sort.bam 
####9. coverage
/home/jupyter/packages/bedtools genomecov -d -ibam example_output_rDNA.bam > example_coverage_test1.txt


