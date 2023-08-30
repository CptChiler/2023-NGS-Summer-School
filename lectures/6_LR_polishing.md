# Workshop: Assembly polishing and variant calling

## Hands-on:

### 1. Polishing

Polish the assembly you produced with `flye`. Use the filtered long-read data. Do this in two steps: using `racon` followed by `medaka`. This will also involve mapping again the long reads to your calculated assemblies! 

#### 1.1. Assembly polishing (Racon)

```bash
# run racon, as input you need the reads, the mapping file, and the assembly you want to polish
racon -t 8 data/ONT_R10-filtered_reads.fastq.gz data/ONT_R10-mapping.sorted.bam flye_output_R10/assembly.fasta > data/ONT_R10-consensus-racon.fasta

# map to new consensus
minimap2 -t 4 -ax map-ont ONT_R10-consensus-racon.fasta data/ONT_R10-filtered_reads.fastq.gz > data/ONT_R10-consensus-racon-mapping.sam
samtools view -bS data/ONT_R10-consensus-racon-mapping.sam | samtools sort -@ 4 > data/ONT_R10-consensus-racon-mapping.sorted.bam  
samtools index data/ONT_R10-racon-onsensus-mapping.sorted.bam

# now look at it in tablet or IGV again. For IGV you have to convert to BAM again and index the mapping file!
igv &

# load assembly file as 'Genomes'
# ->  data/ONT_R10-consensus-racon.fasta

# load mapping file as 'File'
# ->  data/ONT_R10-racon-onsensus-mapping.sorted.bam
```
[Publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411768/) | [Code](https://github.com/isovic/racon)

* a common practice is 2x `racon` polishing followed by 1x `Medaka` (see below)
* for practice, 1x `racon` is fine
* with the newest R10 chemistry and recent basecalling models, it seems that `racon` is also not necessary anymore and people switch to only polish via `Medaka`
* with improvements in ONT accuracy, it might be even the case that long-read polishing is becoming more and more obsolete!

### Assembly polishing and final consensus (Medaka)

`Medaka` is not in your current `workshop` environment because it was conflicting with the other tools. That's why we need a separate Conda environment for `Medaka`:

* Make a new environment for `medaka` 
    * `medaka` might have many dependencies that conflict 
* an alternative to `conda` is `mamba`
    * `mamba` can be much faster in solving your environment, e.g. here for the tool `medaka`
    * thus, let us install `mamba` via `conda` and then install `medaka`

```bash
cd /scratch/$USER/nanopore-workshop
mamba create -y -p envs/medaka "medaka>=1.8.0"
conda activate envs/medaka
```

```bash
# Run Medaka
# ATTENTION: it is always good to assign an appropriate Medaka model -m based on 
# the performed basecalling! Here, we use some example model for a R10.4.1 run with 260 bp/s speed and SUP basecalling. 
# This might to be adjusted based on your data! 
# If you are on the RKI HPC: due to restrictions it might be even difficult to run other Medaka models because 
# they need to be downloaded first.

medaka_consensus -i data/ONT_R10-filtered_reads.fastq.gz -d data/ONT_R10-consensus-racon.fasta \
   -o ONT_R10-medaka -m r1041_e82_260bps_sup_v4.0.0 -t 8

# map to new consensus
minimap2 -t 4 -ax map-ont ONT_R10-medaka/assembly.fasta data/ONT_R10-filtered_reads.fastq.gz > data/ONT_R10-consensus-medaka-mapping.sam
samtools view -bS data/ONT_R10-consensus-medaka-mapping.sam | samtools sort -@ 4 > data/ONT_R10-consensus-medaka-mapping.sorted.bam  
samtools index data/ONT_R10-medaka-onsensus-mapping.sorted.bam

# now look at it in tablet or IGV again. For IGV you have to convert to BAM again and index the mapping file!
igv &

# load assembly file as 'Genomes'
# ->  ONT_R10-medaka/assembly.fasta

# load mapping file as 'File'
# ->  data/ONT_R10-medaka-onsensus-mapping.sorted.bam
```
[Code](https://github.com/nanoporetech/medaka)

**Note** that you should usually change the model parameter (`-m`) to whatever is most appropriate for your basecalling. Also note that `medaka_consensus` is not the same thing as `medaka consensus` (underscore vs space) - the former is a convenience script which does the entire process (including read mapping) while the latter is a subcommand of Medaka which only does the polishing step. (thx to [Ryan Wick for this explanation](https://github.com/rrwick/Trycycler/wiki/Polishing-after-Trycycler)).

Next: [Polishing assembly with short reads](7_SR_polishing.md)
