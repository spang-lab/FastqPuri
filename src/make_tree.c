#include <stdlib.h>
#include <stdio.h>
#include "fa_read.h"

uint64_t alloc_mem; // global variable. Memory allocated in the heap.

int main(){
   char *filename = "/nfs/compdiag/user/pep04706/my_programs/C++/trim_filter/examples/EColi_genome.fa";
   Fa_data *ptr_fa = malloc(sizeof(Fa_data));
   read_fasta(filename,ptr_fa);
//   printf("%s",ptr_fa->entry[0].seq);
   return 0; 
}
