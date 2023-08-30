## Bonus 2

If you also have short-read Illumina data corresponding to your Nanopore data. Use the high-accuracy short-read data to polish your _best_ long-read assembly again. Use `Polypolish` for that. Install `Polypolish` and [familiarize yourself with the tool](https://github.com/rrwick/Polypolish/wiki). For the necessary mapping of the short reads to the assembly you should use `bwa`. Check the `Polypolish` manual!  


```bash
get data
```

### Quick start
(see https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish for more details!)

```sh
bwa index draft.fasta
bwa mem -t 16 -a draft.fasta reads_1.fastq.gz > alignments_1.sam
bwa mem -t 16 -a draft.fasta reads_2.fastq.gz > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish draft.fasta filtered_1.sam filtered_2.sam > polished.fasta
```

Next: [Assembly data analysis](8_Analysis.md)
