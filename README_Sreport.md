# SReport user manual

Creates a summary report in html format of the quality status of a set 
of `fastq` files. Needs the binary output of `Qreport`, for a list 
of files. It will analyze all files with extension `*.bin` inside 
the specified folder in the input and will expect them to 
have the appropriate format (output from `Qreport`).

## Running the program

Usage `C` executable (in folder `bin`): 

```
Usage: ./Sreport -i <INPUT_FOLDER> -o <OUTPUT_FILE> 
Uses all *bin files found in a folder ( output of Qreport) 
and generates a summary report in html format
 -v Prints package version.
 -h Prints help dialog.
 -i Input folder containing *bin data (output from Qreport). Mandatory option.
 -o Output file (with NO extension). Mandatory option.
```


## Output description


- html output:
   * Table with: `# reads , # tiles, % lowQ reads, % reads with N’s` for 
     all samples.
   * Heatmap showing the average quality per position for all samples.

## Example 
 
An example is given in the folder `examples/QReport_Sreport`. To run an 
example, type, 

``` 
    $ cd example/run_test/
    $ ../../../bin/SReport -i ./ -o my_test_summary_report
```
 and compare it with the provided run example, as specified in the README
 file under `./example/QReport_Sreport`

**NOTE:** it has to be run AFTER `Qreport` example

## Contributors

Paula Pérez Rubio 

## License

GPL v3 (see LICENSE.txt)
