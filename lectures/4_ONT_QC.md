# Workshop: Nanopore QC

## Hands-on:

**Note**: If internet connection is slow, we can also distribute the example data via an USB stick (ask your instructor ;) ). 

## 1. Find your data to work

Below are just example folder names, you can also adjust them and use other folder names! Assuming you are on a Linux system on a local machine (laptop, workstation):

```sh
# Switch to a path on your system where you want to store your data and results (you should be already on this path)
cd /scratch/$USER/nanopore-workshop

# make a dir to store your data
mkdir data

# find the path to the stored data
ls /scratch/Tausch/2023-RKI-NGS-Workshop/data

# copy the example data into you working directory
cp -r /scratch/Tausch/2023-RKI-NGS-Workshop/data /scratch/$USER/nanopore-workshop/data

# double-check that everything is in place:
ls -lah data/

# all good? Let's move on to QC!
```

### 1.1 Data managment

After downloading or copying the training data, we will save the path to the respective fastq file in a variable. This is important! We will use from now on this variable to refer to the read file when we start analyzing it (so we can forget about the path). In the code examples the $SAMPLE is the place to put your variable which are listed down below:

```bash
## R10 ONT
raw_ONT_R10 = /scratch/$USER/nanopore-workshop/data/R10.fastq.gz

## Illumina
raw_R1 = /scratch/$USER/nanopore-workshop/data/read_R1.fastq.gz
raw_R2 = /scratch/$USER/nanopore-workshop/data/read_R2.fastq.gz
```

#### Example:
```sh
# To see the first read in R9 fastq file we run:
head -n 4 $raw_ONT_R10

# Note that the $ (dollar sign) makes a word to a variable this wont run:
head -n 4 raw_ONT_R10
```
Always think about how to name your output files. It is advised to keep all infos into your file. For example:
> [Sample]-[what have you done].[extension]
> 
> ONT_R9-filtered_reads.fastq.gz

## 2. ONT basecalled fastq output 

**Attention:** Before doing any calculations it makes sense to start an interactive compute session on the HPC! See [Linux Crash Course](3_setup.md). If you do so, you need to `conda activate envs/workshop` again!

### 2.1. Quality control (NanoPlot)

```bash
cd /scratch/$USER/nanopore-workshop

NanoPlot -t 4 --fastq  $raw_ONT_R10 --title "Raw reads R10" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/raw_R10
```
[Publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty149/4934939) | [Code](https://github.com/wdecoster/NanoPlot)

**Note**: The `\` at the end of a line is only for convenience to write a long command into several lines. It tells the command-line that all lines still belong together although they are separated by "enter" keys. However, if you type all of the command, i.e., paths etc, in one line do not copy/type the backslash at the end of the lines.

### 2.2. Read filtering (Filtlong)

```bash
cd /scratch/$USER/nanopore-workshop

# Note: we use 1 kb as the minimum length cutoff as an example. For your "real" samples other parameters might be better. Do QC before! 
filtlong --min_length 1000 --keep_percent 90 \
    --target_bases 500000000  $raw_ONT_R10 | gzip > data/ONT_R10-filtered_reads.fastq.gz

# Check the quality again:
NanoPlot -t 4 --fastq data/ONT_R10-filtered_reads.fastq.gz --title "Filtered reads R10" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/clean_R10
```
[Code](https://github.com/rrwick/Filtlong)

Next: [Long-read assembly](5_LR_assembly.md)
