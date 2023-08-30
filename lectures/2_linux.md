# Workshop - Linux re-cap

## Short Linux and bash re-cap

* Linux/Bash basics + conda setup ([slides](https://docs.google.com/presentation/d/1Bl-H8lgyPwOUE8EDMSbebaTU14N6QCVTbRyxKrnbTFA/edit?usp=sharing))
* Another good resource: [Introduction to the UNIX command line](https://ngs-docs.github.io/2021-august-remote-computing/introduction-to-the-unix-command-line.html)
* Cheat sheet for Bash: [github.com/RehanSaeed/Bash-Cheat-Sheet](https://github.com/RehanSaeed/Bash-Cheat-Sheet)

A (very short) cheat sheet:
```bash
# Print your user name
echo $USER
# change directory to your user home directory (all of these are the same)
cd /home/$USER
cd $HOME
cd ~  # <- the shortest version, I like this one
# show content of current directory
ls
# make a new directory called 'myfolder'
mkdir myfolder
# make conda environment and activate it
mamba create -n nanoplot
conda activate nanoplot
mamba install nanoplot
# Attention, the above command will create a folder 'nanoplot' in your default path, e.g. ~/miniconda3/envs
# However, you can also specify any other folder, and which we will also do in the training:
mamba create -p envs/nanoplot
conda activate envs/nanoplot
mamba install nanoplot
# run a program
NanoPlot reads.fq.gz ...
```

Next: [Linux for Bioinformatics re-cap](lectures/3_setup.md)