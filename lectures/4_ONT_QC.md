# Workshop: Nanopore QC

**Note**: If internet connection is slow, we can also distribute the example data via an USB stick (ask your instructor ;) ). 

## 1. Create a folder for the hands-on work

Below are just example folder names, you can also adjust them and use other folder names! Assuming you are on a Linux system on a local machine (laptop, workstation):

```sh
# Switch to a path on your system where you want to store your data and results (you should be already on this path)
cd /scratch/$USER
# Create your workspace folder you want to work
mkdir nanopore-workshop
cd nanopore-workshop
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

Nanopore data can be basecalled during the run. You can also run basecalling yourself, but it is time consuming and we'll skip it here. The output is sorted into barcodes and in every folder are a list of compressed fastq files. These files need to be combined into one file for further use.

```bash
# make sure that you are located in your workshop folder, if not cd
cat <run name>/fastq_pass/*.fastq.gz > <file name>.fastq.gz
```

After downloading or copying the training data, we will save the path to the respective fastq file in a variable called SAMPLE. This is important! We will use from now on this variable to refer to the read file when we start analyzing it (so we can forget about the path).

```bash
## R9 ONT
raw_ONT_R9 =

## R10 ONT
raw_ONT_R10 = 

## Illumina
raw_R1 =
raw_R2 =
```

## 2. ONT basecalled fastq output 

**Attention:** Before doing any calculations it makes sense to start an interactive compute session on the HPC! See [Linux Crash Course](3_setup.md). If you do so, you need to `conda activate envs/workshop` again!

### 2.1. Quality control (NanoPlot)

```bash
cd /scratch/$USER/nanopore-workshop
NanoPlot -t 4 --fastq data/barcode01.fastq.gz --title "Raw reads" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/raw
```
[Publication](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty149/4934939) | [Code](https://github.com/wdecoster/NanoPlot)

**Note**: The `\` at the end of a line is only for convenience to write a long command into several lines. It tells the command-line that all lines still belong together although they are separated by "enter" keys. However, if you type all of the command, i.e., paths etc, in one line do not copy/type the backslash at the end of the lines.

### 2.2. Read filtering (Filtlong)

```bash
# Note: we use 1 kb as the minimum length cutoff as an example. For your "real" samples other parameters might be better. Do QC before! 
filtlong --min_length 1000 --keep_percent 90 \
    --target_bases 500000000 data/barcode01.fastq.gz > barcode01-filtered.fastq

# Check the quality again:
NanoPlot -t 4 --fastq barcode01-filtered.fastq --title "Filtered reads" \
    --color darkslategrey --N50 --loglength -f png -o nanoplot/clean
```
[Code](https://github.com/rrwick/Filtlong)

Next: [Long-read assembly](5_LR_assembly.md)
