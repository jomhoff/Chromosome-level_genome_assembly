
## Using Hi-C seqeuncing to create a chromosome-level genome 

For this genome, I started with my long-read hifiasm assembly. The assembly methods can be found [here](https://github.com/jomhoff/Genome-Assembly).
Hi-C sequencing was done through Phase Genomics.

## Align Hi-C data to Assembly with BWA

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

## Remove duplicates and sort 

```
#!/bin/sh
#SBATCH --job-name=stools
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=25
#SBATCH --mem=50gb
#SBATCH --time=100:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=jhoffman1@amnh.org

conda activate hic
cd /home/jhoffman1/mendel-nas1/Hi-C

samtools sort -@ 20 -n -o hic_reads.querysorted.bam hic_reads.bam &
wait
samtools fixmate -m hic_reads.querysorted.bam hic_reads.fixmate.bam &
wait
samtools markdup -@ 20 -r hic_reads.fixmate.sorted.bam hic_reads.dedup.bam &
wait
samtools index -@ 20 hic_reads.dedup.bam
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

./yahs /home/jhoffman1/mendel-nas1/Hi-C/juicer/references/plestiodonFasciatus.softmasked_sf.fasta /home/jhoffman1/mendel-nas1/Hi-C/hic_reads.dedup.bam
```

## Generate Hi-C contact map 

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

```
#!/bin/bash
#SBATCH --job-name=juicerpre
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=50G
#SBATCH --time=48:00:00
#SBATCH --output=juicerpre_%j.out
#SBATCH --error=juicerpre_%j.err

conda activate /mendel-nas1/jhoffman1/miniforge3/envs/hic
module load Java/jdk-1.8.0_281

java -Xmx50G -jar juicer_tools_1.22.01.jar pre trimmed_file.txt pfas.hic genome.chrom.sizes -v -j 10
```
