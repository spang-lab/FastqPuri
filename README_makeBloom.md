# makeTree user manual

Reads a `fasta` file, creates a bloom filter with a predefined:
 - size (bits). 
 - number of hash functions.
 - false positive rate.
and saves it to a file. 

## Running the program

Usage `C` executable (in folder `bin`): 

```
Usage: ./makeBloom --fasta <FASTA_INPUT> --output <FILTERFILE> --kmersize [KMERSIZE] 
 (--fal_pos_rate [p] | --hashNum [HASHNUM] | --bfsizeBits [SIZEBITS])
Options: 
 -v, --version      Prints package version.
 -h, --help         Prints help dialog.
 -f, --fasta        Fasta input file. Mandatory option.
 -o, --output       Output file, with NO extension. Mandatory option.
 -k, --kmersize     kmer size, number or elements. Optional(default = 25)
 -p, --fal_pos_rate false positive rate. Optional (default = 0.05)
 -g, --hashNum      number of hash functions used. Optional (default
                    value computed from the false positive rate).
 -m, --bfsizeBits   size of the filter in bits. It will be forced to be
                    a multiple of 8. Optional (default value computed
                    from the false positive rate).
NOTE: the options -p, -n, -m are mutually exclusive. The program 
      will give an error if more than one of them are passed as input.
      It is recommended to pass the false positive rate and let the 
      program compute the other variables (excepting singular situations)
      If none of them are passed, the false positive rate will default
      to 0.05 and the other variables will be computed respecting this
      requirement. See documentation and supplementary material for 
      more details.
```

## Generating a bloom filter

TODO

### Input parameters. 

### Encoding the kmers

### Incluing a kmer in the filter

For further details, read the `Doxygen` documentation of the files
`bloom.c`, `init_makeBloom.c`, `makeBloom.c` and the supplementary material.


## Output description

Two files are created as output: 
 - `<FILTERFILE>.bf`: contains the filter (binary file).
 - `<FILTERFILE>.txt`: contains the following information: 
   * kmersize: length of the kmers inserted in the filter, 
   * bfsizeBits: size of Bloom filter in bits,
   * hashNum: number of hash functions used in the filter, 
   * falsePosRate: false positive rate,
   * nelem number: of elemens (kmers in the sequece) contained in the filter.

For further details, read the `Doxygen` documentation of the files
`bloom.c`, `init_makeBloom.c`, `makeBloom.c`

## Example 
 
TODO

## Contributors

Paula Pérez Rubio 

## License

GPL v3 (see LICENSE.txt)
