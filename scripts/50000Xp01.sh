#!/bin/bash
#SBATCH --job-name=50000Xp01
#SBATCH --output=50000Xp01.out
#SBATCH --error=50000Xp01.err
#SBATCH --verbose
#SBATCH --array=1
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=36GB
#SBATCH --mail-type=END
#SBATCH --mail-user=xm595@nyu.edu

#module purge
#module load
module purge
#module load samtools/intel/1.14
#variant frequency = 0.01
#control parameters: c (50000*0.01 == 500), M, and o

module load bowtie2/2.4.4
module load lofreq/2.1.5

#this module swap is important to run before loading samtools(suspect it has something to do with the multiple versions available
module swap htslib/intel/1.12 htslib/intel/1.14

module load samtools/intel/1.14
module load bedtools/intel/2.29.2

#singularity exec --nv \
#            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
#            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
#            /bin/bash -c "source /ext3/env.sh;conda activate py27; python neat-genreads/utilities/computeGC.py -r /scratch/cgsb/hochwagen/Human_rDNA_project/syn2/rDNA_prototype_prerRNA_only.fa -i /scratch/cgsb/hochwagen/Human_rDNA_project/synetic_data/ERR3240249_rDNA.coverage -o gcmodel.p"

#samtools view /scratch/cgsb/hochwagen/Human_rDNA_project/synetic_data/ERR3240249_sort.bam | python neat-genreads-master/utilities/computeFraglen.py

#singularity exec --nv \
#            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
#            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
#            /bin/bash -c "source /ext3/env.sh;conda activate py27; module load samtools/intel/1.14; samtools view /scratch/cgsb/hochwagen/Human_rDNA_project/synetic_data/ERR3240249_sort.bam | python neat-genreads-master/utilities/computeFraglen.py"



#singularity exec --nv \
#            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
#            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
#            /bin/bash -c "source /ext3/env.sh;conda activate py27; python neat-genreads/utilities/genSeqErrorModel.py -i /scratch/cgsb/hochwagen/Human_rDNA_project/synetic_data/ERR3240249_1.fastq -i2 /scratch/cgsb/hochwagen/Human_rDNA_project/synetic_data/ERR3240249_2.fastq -o seq_error.p" 

#######start

singularity exec --nv \
            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
            /bin/bash -c "source /ext3/env.sh;conda activate py27; python neat-genreads/genReads.py -r /scratch/cgsb/hochwagen/Human_rDNA_project/syn2/rDNA_prototype_prerRNA_only.fa -R 150 -o out_50000Xp01_1 --bam --vcf --pe-model fraglen.p -e seq_error.p --gc-model gcmodel.p -p 1 -M 0.015 -c 500 -t human_rDNA_repeat_benchmark.bed -to 0.4 --rng 123"

singularity exec --nv \
            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
            /bin/bash -c "source /ext3/env.sh;conda activate py27; python neat-genreads/genReads.py -r /scratch/cgsb/hochwagen/Human_rDNA_project/syn2/rDNA_prototype_prerRNA_only.fa -R 150 -o out_50000Xp01_2 --bam --pe-model fraglen.p -e seq_error.p --gc-model gcmodel.p -p 1 -M 0 -c 49500 -t human_rDNA_repeat_benchmark.bed -to 0.4 --rng 456"

singularity exec --nv \
            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
            /bin/bash -c "source /ext3/env.sh;conda activate py27;module load samtools/intel/1.14; python neat-genreads/mergeJobs.py -i out_50000Xp01_1 out_50000Xp01_2 -o simulation_50000Xp01 -s /share/apps/samtools/1.14/intel/bin/samtools"

sample_list=("simulation_50000Xp01_read")
for sample in "${sample_list[@]}"
do
    start1=$(date +%s.%N)
    echo "Processing sample: $sample"

    #fetchNdump
    
#    prefetch ${sample}
#    fasterq-dump ${sample}
#    rm -r ${sample}

    #align reads to rDNA
    bowtie2 -5 1 -N 1 -p 8 \
    -x /scratch/cgsb/hochwagen/Human_rDNA_project/rDNA_prototype_prerRNA_only/rDNA_prototype_prerRNA_only \
    -1 ${sample}1.fq \
    -2 ${sample}2.fq \
    -S ${sample}_output.sam
#    rm ${sample}_1.fastq
#    rm ${sample}_2.fastq

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
#    rm ${sample}_rDNA.bam
####move files to it corresponding folder.
    
    mkdir H_${sample}
#    mv ${sample}_rDNA.vcf ${sample}/
    mv ${sample}_rDNA.bam H_${sample}/
    mv ${sample}_rDNA.cram H_${sample}/
    mv ${sample}_rDNA_coverage.txt H_${sample}/
    echo "Processing of sample $sample complete"
    dur1=$(echo "$(date +%s.%N) - $start1" | bc)
    printf "Execution time: %.6f seconds" $dur1
done



singularity exec --nv \
            --overlay /scratch/xm595/pytorch-example/overlay-15GB-500K.ext3:rw \
            /scratch/work/public/singularity/cuda11.6.124-cudnn8.4.0.27-devel-ubuntu20.04.4.sif\
            /bin/bash -c "source /ext3/env.sh;conda activate py27;python neat-genreads/utilities/vcf_compare_OLD.py -r /scratch/cgsb/hochwagen/Human_rDNA_project/syn2/rDNA_prototype_prerRNA_only.fa -g simulation_50000Xp01_golden.vcf -w simulation_50000Xp01_read_rDNA.vcf -o simulation_50000Xp01 -a 0.002 --vcf-out --incl-fail --no-plot"

mv out_50000Xp01_* H_simulation_50000Xp01_read
mv simulation_50000Xp01* H_simulation_50000Xp01_read
