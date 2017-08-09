FastqArazketa, a fq quality control and filter tool 
=========

Software and source code of `FastqArazketa`. Creates quality 
reports of `fastq` files and filters them removing low quality 
reads, reads containing too many N's or contamination reads 
(unwanted rRNA reads, impurities coming from another organism, ...).


## Installation

Clone the repository, or download the source. Make sure that 
your system supplies the following dependencies for FastqArazketa.


- cmake 
- a C compiler supporting the c11 standard (change the compiler flags otherwise)
- pandoc (optional, see documentation in PANDOC.md)
- Rscript (optional)
- Following R packages installed (optional):
   * pheatmap
   * knitr
   * rmarkdown

**NOTE:**  FastqArazketa will work without the optional dependencies 
but will skip creating html reports if they are not available.

```
cmake -H. -Bbuild/ [-DRSCRIPT=myRscriptexec ... ]
cd build 
make 
```

When running cmake, there are some variables you can set 
using the option -D followed by the variable name. This variables are, 

- CMAKE_C_COMPILER: C compiler (default gcc)
- CMAKE_C_FLAGS: compiler flags (default `-Wall -O3 -march=native -std=c11`).
- PANDOC: pandoc executable (default `pandoc`),
- RSCRIPT:  Rscript executable (default `Rscript`), 
- READ_MAXLEN: Maximum Illumina read length (default 400),
- RMD_QUALITY_REPORT path to `quality_report.Rmd`,
- RMD_SUMMARY_REPORT path to `summary_report.Rmd`,

The executables will be created in the folder `bin`. 

# Executables

* `Qreport`: generates a quality report in html format (see README_Qreport.md),
* `Sreport`: generates a summary report in html format (see README_Sreport.md),
* `makeSA`: work in progress,
* `makeBloomFilter`: work in progress,
* `makeTree`: work in progress,
* `trim_filter`: work in progress.

# Documentation of the code

A Doxygen documentation of the code is available: 
- `html` version under the folder `html` (open `index.html` with a browser).
- `pdf` version: `latex/refman.pdf`
