FastqPuri, an fq quality control and filter tool 
=========

Software and source code of `FastqPuri`. It creates quality reports of
`fastq` files and filters them removing low quality reads, reads
containing too many N's or contamination reads (unwanted rRNA reads,
impurities coming from another organism, ...).


## Installation

Clone the repository, or download the source. Make sure that 
your system supplies the following dependencies for FastqPuri.

- OS: Linux (clang, gcc), Mac OS (clang, gcc), OpenBSD (clang)
- `cmake` (at least version 2.8), 
- a `C` compiler supporting the `c11` standard 
  (change the compiler flags otherwise),
- pandoc (optional, see documentation in `PANDOC.md`),
- `Rscript` (optional),
- Following `R` packages installed (optional):
   * `pheatmap`
   * `knitr`
   * `rmarkdown`

**NOTE:**  FastqPuri will work without the optional dependencies 
but will skip creating html reports if they are not available.

```
$ cmake -H. -Bbuild/ [-DRSCRIPT=/path/to/my/R/bin/Rscript] [-DCMAKE_INSTALL_PREFIX=/path/to/my/root] ... 
$ cd build 
$ make 
$ sudo make install  
```

When running `cmake`, there are some variables you can set 
using the option -D followed by the variable name. These variables are:

- `CMAKE_C_COMPILER`: `C` compiler (default `gcc`)
- `CMAKE_C_FLAGS`: compiler flags (default `-Wall -O3 -march=native -std=c11`).
- `CMAKE_INSTALL_PREFIX`: root path for `make install`, e.g. to
  redirect to a directory with user access (default /usr/local),
- `PANDOC`: `pandoc` executable (default `pandoc`),
- `RSCRIPT`: `Rscript` executable (default `Rscript`),
- `READ_MAXLEN`: Maximum Illumina read length
- (default 400),

The executables will be created in the folder `bin` and installed in `/usr/local/bin`. 
`R` scripts will be installed in `/usr/local/share/FastqPuri/R`. 

**WARNING:** do not move the executables that depend on `R` scripts, 
anywhere else, unless you also move the corresponding `R` scripts respecting
the local folder structure. 


## Executables

* `Qreport`: creates a quality report in html format (see `README_Qreport.md`),
* `Sreport`: creates a summary report in html format on a set of samples, 
   regarding either the original files or the filtering process
   (see `README_Sreport.md`),
* `makeBloom`: creates a  bloom filter from a fasta file of a certain size,
   and stores it in a file (see `README_makeBloom.md`)
* `makeTree`: creates a tree of a certain depth from a fasta file and stores
   it in a file (see `README_makeTree.md`),
* `trimFilter`: performs the filtering process for single-end data 
   (see `README_trimFilter.md`).
* `trimFilterPE`: performs the filtering process for double stranded data 
   (see `README_trimFilterPE.md`).

An exemplar work flow could be:

* `Qreport`
* `Sreport`
* `makeBloom`
* `trimFilter` or `trimFilterPE`
* `Qreport`
* `Sreport`

## Documentation of the code

A Doxygen documentation of the code is available:

- `html` version under the folder `html` (open `index.html` with a browser).
- `pdf` version: `latex/refman.pdf`

## Use a docker container to run FastqPuri

The file 'Dockerfile' documents the exact linux installation we used
for testing. If you have a docker installation ready on your machine,
you may want to use a docker container for easy installation and
capsulated usage of FastqPuri. After cloning this project from github
and change to its main directory, you may install a docker container
as follows:

```
$ docker build -t fastqpuri .
```

This will create a container based on the debian linux distribution
covering all dependencies including R and pandoc.  As soon as such a
container is installed, you can use it either interactively:

```
$ docker run -v $PWD:/tmp -it fastqpuri
```

or by running a pipeline implemented in an executable bash script:

```
$ docker run -v $PWD:/tmp fastqpuri ./pipeline.sh
```

Note that this call generates results in the docker container
directory `/tmp` but also keeps them after closing the docker container
locally where the container was started.

Instead of generating the docker container yourself with 'docker
build', you can also pull a pre-built image from the docker hub as
follows:

```
$ docker pull clottaz/fastqpuri
```

You can run such a pre-built image with 'docker run' by indicating the
images as 'clottaz/fastqpuri'.


## Use a singularity container to run FastqPuri

Alternativly, if you have singularity installed on your machine, you
can call our docker container for FastqPuri as follows:

```
$ singularity shell --bind .:/tmp docker://clottaz/fastqpuri
```

This call opens a shell within the container.
With `--bind` we  mount the current directory also in the container.
The syntax is as follows: --bind src:dest; src is the source path on
the host and dest is the destination path in the container, i.e. where
you would like to make the source path available in your container.
Note that this destination path in your container should be an existing
directory, the operation will fail if you do not create the directory first.
Hence, when we call `singularity shell` like this, the working directory
in the container is `/tmp`.

Alternatively, in order to execute a script from the current
directory, call singularity as follows:

```
$ singularity run --bind .:/tmp docker://clottaz/fastqpuri /tmp/pipeline.sh
```

Note that `/tmp/pipeline.sh` relates to the call within the
container. Thus, `pipeline.sh` is located in the directory where singularity
run is executed, but will be made available to the container via the `--bind`
parameter.

If you want to invoke a function of FastqPuri, you can use the 'exec'
command like so:

```
singularity exec docker://clottaz/fastqpuri Qreport -h
```

or invoke a script located in your home directory (assuming that
run_ex_TREE.sh is located in your home directory):

```
$ singularity exec docker://clottaz/fastqpuri $HOME/run_ex_TREE.sh
```

Singularity documentation can be found here: [https://www.sylabs.io/docs/](https://www.sylabs.io/docs/)



## Installation via bioconda **-under construction**.

*We are currently working on a bioconda environment for FastqPuri. 
If you follow the instructions below, it is quite likely that 
FastqPuri will not yet properly run from the bioconda environment.
Sorry about that and please stay tuned!*

Bioconda is a channel for the conda package manager specializing in
bioinformatics software. Have a look at the reference:

* Bjoern Gruening, Ryan Dale, Andreas Sjoedin, Brad A. Chapman, Jillian
  Rowe, Christopher H. Tomkins-Tinch, Renan Valieris, the Bioconda
  Team, and Johannes Koester. 2018. Bioconda: Sustainable and
  Comprehensive Software Distribution for the Life Sciences. Nature
  Methods, 2018.

To find out how to use bioconda, see [https://bioconda.github.io](https://bioconda.github.io).
For installing FastqPuri in a bioconda environment, you have to install
either `miniconda` or `anaconda` and register channels as follows:

```
$ conda config --add channels defaults
$ conda config --add channels bioconda
$ conda config --add channels conda-forge
```

Then you can install `fastqpuri`:

```
$ conda install fastqpuri
```

Actually, you may also want to use a specific environment for the
sequencing quality control:

```
$ conda create -n qc fastqpuri
```

This call installs `FastqPuri` directly in a separate environment.

## Contributors

Paula Pérez Rubio,
Claudio Lottaz,
Julia Engelmann 

## License

GPL v3 (see LICENSE.txt)
