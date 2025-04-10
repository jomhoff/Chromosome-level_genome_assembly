
## Using Hi-C seqeuncing to create a chromosome-level genome 

For this genome, I started with my long-read hifiasm assembly. The assembly methods can be found [here](https://github.com/jomhoff/Genome-Assembly).
Hi-C sequencing was done through Phase Genomics.

## Align Hi-C data to Assembly with BWA
make sure to index draft genome with BWA before starting
```
bwa index hoff_pfas_total.fasta
```
Then run BWA
```
#!/bin/sh
#SBATCH --job-name=bwa
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=50gb
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org

conda activate /mendel-nas1/jhoffman1/miniforge3/envs/hic
cd /home/jhoffman1/mendel-nas1/Hi-C/nonmasked

DRAFT="/home/jhoffman1/mendel-nas1/fasciatus_genome/KY_LR/assembly/hoff_pfas_total.fasta"
HICR1="/home/jhoffman1/mendel-nas1/Hi-C/juicer/work/pfas/fastq/amnh-fs20603_HiC_R1.fastq.gz"
HICR2="/home/jhoffman1/mendel-nas1/Hi-C/juicer/work/pfas/fastq/amnh-fs20603_HiC_R2.fastq.gz"

bwa mem -t 25 -5SP $DRAFT $HICR1 $HICR2 | samtools view -bS - > hic_reads.bam
```

## Remove duplicates and sort 

```
#!/bin/sh
#SBATCH --job-name=stools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=50gb
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org

conda init
conda activate /mendel-nas1/jhoffman1/miniforge3/envs/hic
cd /home/jhoffman1/mendel-nas1/Hi-C
# Name sort
samtools sort -n -@ 20 -o hic_reads.qnamesorted.bam hic_reads.bam
wait
# Fixmate
samtools fixmate -m -@ 20 hic_reads.qnamesorted.bam hic_reads.fixmate.bam
wait
# Coordinate sort 
samtools sort -@ 20 -o hic_reads.fixmate.sorted.bam hic_reads.fixmate.bam
wait
# Index 
samtools index -@ 20 hic_reads.fixmate.sorted.bam
wait
# Remove duplicates
samtools markdup -@ 20 -r hic_reads.fixmate.sorted.bam hic_reads.dedup.bam
wait
# Name sort for YAHS
samtools sort -n -@ 20 -o hic_reads.final.bam hic_reads.dedup.bam
```


## Scaffolding with [YAHS](https://github.com/c-zhou/yahs)

```
#!/bin/bash
#SBATCH --job-name=yahs
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --time=10:00:00
#SBATCH --output=yahs_%j.out
#SBATCH --error=yahs_%j.err

conda activate /mendel-nas1/jhoffman1/miniforge3/envs/hic
cd /home/jhoffman1/mendel-nas1/Hi-C/yahs

./yahs /home/jhoffman1/mendel-nas1/fasciatus_genome/KY_LR/assembly/hoff_pfas_total.fasta /home/jhoffman1/mendel-nas1/Hi-C/hic_reads.final.bam -o ./nonmasked/pfas.out
```

## Generate Hi-C contact map 
Before starting, I renamed the final fasta scaffolds_final.fa
Also, I indexed my draft genome fasta with
```
samtools faidx hoff_pfas_total.fasta
```

Now, create the input .txt file for juicer from yahs output (provided by yahs)
```
#!/bin/bash
#SBATCH --job-name=yahsj
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=juc_%j.out
#SBATCH --error=juc_%j.err

/home/jhoffman1/mendel-nas1/Hi-C/yahs/juicer pre /home/jhoffman1/mendel-nas1/Hi-C/hic_reads.final.bam pfas.out_scaffolds_final.agp /home/jhoffman1/mendel-nas1/fasciatus_genome/KY_LR/assembly/hoff_pfas_total.fasta.fai \
| sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G \
| awk 'NF' > alignments_sorted.txt.part && mv alignments_sorted.txt.part alignments_sorted.txt
```

Used this code to generate the required chromosome size file used to make a Hi-C contact map (.hic)
```
samtools faidx pfas.out_scaffolds_final.fa
cut -f1,2 pfas.out_scaffolds_final.fa.fai > scaffolds_final.chrom.sizes
```

Now, use juicertools to generate the .hic file
```
#!/bin/bash
#SBATCH --job-name=contactmap
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=50G
#SBATCH --time=10:00:00
#SBATCH --output=cmap_%j.out
#SBATCH --error=cmap_%j.err

conda activate /mendel-nas1/jhoffman1/miniforge3/envs/hic
module load Java/jdk-1.8.0_281

cd /home/jhoffman1/mendel-nas1/Hi-C/yahs/nonmasked

java -Xmx50G -jar /home/jhoffman1/mendel-nas1/Hi-C/yahs/juicer_tools_1.22.01.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && 
mv out.hic.part nonmasked_out.hic
```
The resulting .hic can be visualized with Juicebox Assembly Tools
