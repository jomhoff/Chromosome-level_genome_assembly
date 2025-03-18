
## Using Hi-C seqeuncing to create a chromosome-level genome (following the general pipeline from the [Genome Assembly Cookbook](https://aidenlab.org/assembly/manual_180322.pdf))

For this genome, I started with my long-read hifiasm assembly. The assembly methods can be found [here](https://github.com/jomhoff/Genome-Assembly).
Hi-C sequencing was done through Phase Genomics.

## Align Hi-C data to Assembly with BWA and sort the resulting .bam file

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
amnh-fs20603_1504586_S3HiC_R1.fastq.gz \
amnh-fs20603_1504586_S3HiC_R2.fastq.gz | \
samtools sort -@ 25 -o hic_reads.sorted.bam
```

