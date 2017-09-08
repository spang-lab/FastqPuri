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

## Details on bloom filters

**NOTE**: For further details, read the `Doxygen` documentation of the files
`bloom.c`, `init_makeBloom.c`, `makeBloom.c` and the supplementary material.

A bloom filter is a probabilistic data structure used to test if an element
is a member of a set. False positive matches are possible but false negative 
are not. For a given set of `n` elements, we proceed as follows: 

- decide on the number `g` of hash functions we will use for the construction
of the filter and the length of the filter `m` (number of bits). This choice 
will be made based on the false positive rate we want to 
achieve (see Parameters)

- we create an empty bloom filter, `B`, i.e. and array of `m` bits set 
to `0`. For every element in the set,  **s<sub>&alpha;</sub> &isin; S**, the 
`g` hash functions,  **H<sub>i</sub> (s<sub>&alpha;</sub>) &forall; 
i  &isin; {1,...,g}** and set the corresponding bits to
`1` in the filter, i.e., 
**B[H<sub>i</sub> (s<sub>&alpha;</sub>) mod m] = 1 &forall; i &isin; 
{1,...,g}** and `0` otherwise. 

Then, if we want to check whether an element **s<sub>&beta;</sub>** is in 
the set **S**, we compute **H<sub>i</sub> (s<sub>&beta;</sub>) &forall;
i &isin; {1,...,g}** and check whether all coresponding positions in the 
filter are set to `1`, in which case we can say that **s<sub>&beta;</sub>** 
might be in the set. Otherwise it is definitely not in the set. 

### Parameters. 

We choose the parameters so that the desired false positive rate is 
achieved. Alternatively, we can pass the filter size, and then the 
number of hash functions to be used is tuned so that the false positive
rate is minimized. 

We assume the hash functions select all positions with the same 
probability. The probability that an bit in the filter `B` is not 
set to `1` after inserting an element using `g` hash functions is: 

<p align="center">
**(1 - <sup>1</sup>&frasl;<sub>m</sub>)<sup>g</sup>**
</p>

where `m` is the number of bits of the filter. If we insert `n` elements, 
the probability that an element is still 0 is: 

<p align="center">
**(1 - <sup>1</sup>&frasl;<sub>m</sub>)<sup>gn</sup>**
</p>

The probability that a bit is 1 is then, 


### Creating a bloom filter from a fasta file

### Checking ir a read in a fastq file is in the filter



## Contributors

Paula Pérez Rubio 

## License

GPL v3 (see LICENSE.txt)
