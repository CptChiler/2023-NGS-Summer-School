## (H)PC setup

If you are using your own (Linux) laptop or one provided for the workshops you are good to go. If you are using an RKI Windows laptop, you have to connect to the Linux High-Performance Cluster (HPC). You also need an account that you can request via an IT ticket. 

### HPC access

* Install `MobaXterm` via the RKI Software Kiosk
* Open `MobaXterm`
* Connect to the HPC login node, your instructors will tell you the name
    * Select "Session": "SSH"
    * "Remote host": "provided login node name" 
    * "Username": RKI account name
 
### Filezilla access

* Install `Filezilla` via the RKI Software Kiosk
* Detailed information can be found in the RKI Confluence, search for:"
    * "HPC-Datentransfer mit FileZilla"
  
### HPC usage

* Detailed information on HPC infrastructure and usage can be found in the RKI Confluence, search for:
    * "HPC Aufbau"
    * "HPC Nutzung"
    * "HPC FAQ"
* Opening an interactive Shell
    * On the HPC we have login and compute nodes. We dont want to compute on login nodes.
    * An interactive shell is simply any shell process that you use to type commands, and get back output from those commands. That is, a shell with which you interact. We want to connect to a compute node for the workshop.
    * Open `MobaXterm`` and connect to one of the login nodes (ask instructors)

Opening an interactive shell on the RKI HPC:
```sh
#start an interactive bash session using the default ressources
srun --pty bash -i

#start an interactive bash session using 8 CPUs, 32GB RAM, 30GB HDD reserved for 8 hours
srun --time=6:00:00 --cpus-per-task=8 --mem=32GB --gres=local:30 --pty bash -i

#IMPORTANT to free the blocked resources after our work is done close the interactive shell type:
exit
```
Due to competing requests it may take some time until the requested resources can be provided by the system. Therefore, wait patiently until the prompt appears. Reducing requested resources might help as well.

## Prepare working directory and prepare mamba/conda 

* Mamba is a packaging manager that will help us to install bioinformatics tools and to handle their dependencies automatically
* Mamba works together with the conda package manager, and makes installing packages faster
* You will use the mamba command to create environments and install packages, and conda command for some other package management tasks like configuration and activating environments (yes it can be a bit confusing)
* In the terminal enter:

```bash
# Switch to a directory with enough space
cd /scratch/$USER

# make a new folder called 'nanopore-workshop' (this is where we will work in)
mkdir nanopore-workshop

# switch to this folder
cd  nanopore-workshop
```

* Set up mamba

```bash
# Test if conda and mamba is installed (should be rdy on HPC)
conda --help
mamba --help

# ATTENTION: the space in your home directory might be limited (e.g. 10 GB) and per default conda installs tools into ~/.conda/envs
# Thus, take care of your disk space!
# On the HPC you can take care of this by moving ~/.conda to /scratch and making a symlink from your home directory:
# mv ~/.conda /scratch/dot-conda
# ln -s /scratch/dot-conda ~/.conda

# You should now be able to create environments, install tools and run them

# add repository channels for bioconda
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
```

* Create and activate a new conda environment

```bash
# Working directory
cd /scratch/$USER/nanopore-workshop

# create folder for enviroments
# We want to save all programms locally to save space
mkdir envs

# -n parameter to specify the name
mamba create -p envs/workshop

# activate this environment
conda activate envs/workshop

# You should now see (workshop) at the start of each line.
# You switched from the default 'base' environment to the 'workshop' environment.
# Which is placed in a folder envs/workshop
```

Next: [Long-read Nanopore Introduction & Quality Control](4_ONT_QC.md)
