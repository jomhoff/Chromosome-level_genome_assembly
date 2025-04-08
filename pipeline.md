
## Using Hi-C seqeuncing to create a chromosome-level genome 

For this genome, I started with my long-read hifiasm assembly. The assembly methods can be found [here](https://github.com/jomhoff/Genome-Assembly).
Hi-C sequencing was done through Phase Genomics.

## Align Hi-C data to Assembly with BWA
make sure to index draft genome with BWA before starting
```
bwa index hoff_pfas_total.fasta
```
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

conda activate hic
cd /home/jhoffman1/mendel-nas1/Hi-C

bwa mem -t 25 -5SP plestiodonFasciatus.softmasked_sf.fasta \
amnh-fs20603_HiC_R1.fastq.gz \
amnh-fs20603_HiC_R2.fastq.gz | samtools view -bS - > hic_reads.bam
```
or 

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

conda activate hic
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
#SBATCH --mem=2gb
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org

conda activate hic
cd /home/jhoffman1/mendel-nas1/Hi-C

samtools sort -@ 20 -n -o hic_reads.querysorted.bam hic_reads.bam &
wait
samtools index -@ 20 hic_reads.dedup.bam hic_reads.querysorted.bam
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

./yahs /home/jhoffman1/mendel-nas1/fasciatus_genome/KY_LR/assembly/hoff_pfas_total.fasta /home/jhoffman1/mendel-nas1/Hi-C/hic_reads.dedup.bam
```

## Generate Hi-C contact map 
Before starting, I renamed the final fasta scaffolds_final.fa
Also, I indexed my draft genome fasta with
```
samtools faidx plestiodonFasciatus.softmasked_sf.fasta
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

./juicer pre hic-to-contigs.bin scaffolds_final.agp /home/jhoffman1/mendel-nas1/Hi-C/juicer/references/plestiodonFasciatus.softmasked_sf.fasta.fai \
| sort -k2,2d -k6,6d -T ./ --parallel=8 -S32G \
| awk 'NF' > alignments_sorted.txt.part && mv alignments_sorted.txt.part alignments_sorted.txt
```

Used this code to generate the required chromosome size file used to make a Hi-C contact map (.hic)
```
samtools faidx scaffolds_final.fa
cut -f1,2 scaffolds_final.fa.fai > scaffolds_final.chrom.sizes
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

cd /home/jhoffman1/mendel-nas1/Hi-C/yahs

java -Xmx50G -jar juicer_tools_1.22.01.jar pre alignments_sorted.txt out.hic.part scaffolds_final.chrom.sizes && 
mv out.hic.part out.hic
```
