/****************************************************************************
 * Copyright (C) 2017 by Paula Perez Rubio                                  *
 *                                                                          *
 * This file is part of FastqArazketa.                                      *
 *                                                                          *
 *   FastqArazketa is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as                *
 *   published by the Free Software Foundation, either version 3 of the     *
 *   License, or (at your option) any later version.                        *
 *                                                                          *
 *   FastqArazketa is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with FastqArazketa.                                              *
 *   If not, see <http://www.gnu.org/licenses/>.                            *
 ****************************************************************************/

/**
 * @file trimFilterDS.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 25.08.2017
 * @brief trimFilterDS main function
 *
 * This file contains the trimFilterDS main function.
 * See README_trimFilterDS.md and README_trimFilter for more details.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "defines.h"
#include "config.h"
#include "fopen_gen.h" 
#include "ds_read.h"
#include "fq_read.h"
#include "init_trimFilter.h" // change in the future but ok for our current tests

uint64_t alloc_mem = 0;  /**< global variable. Memory allocated in the heap.*/
Iparam_trimFilter par_TF;  /**< global variable: Input parameters of makeTree.*/

int main(int argc, char *argv[]) {
  char fq_r1[MAX_FILENAME] = "../../../../programming_tests/C++_tests/aligner/human_reads.bwa.read1.fastq.gz";
  char fq_r2[MAX_FILENAME] = "../../../../programming_tests/C++_tests/aligner/human_reads.bwa.read2.fastq.gz"; 
  par_TF.L = 50; 
  FILE *f_r1, *f_r2; 
  f_r1 = fopen_gen(fq_r1, "r");
  f_r2 = fopen_gen(fq_r2, "r");
  int newl1=0, newl2=0; 
  int offset1 = 0, offset2 = 0;
  int l1_i=0, l1_f=0, l2_i=0, l2_f=0; 
  int j1=0, j2=0; 
  int nl1=0, nl2=0; 
  int stop1=0, stop2=0;  
  int nentries = 0; 
  char *buffer1 = malloc(sizeof(char)*(B_LEN + 1));
  char *buffer2 = malloc(sizeof(char)*(B_LEN + 1));
  Fq_read  *seq1 = malloc(sizeof(Fq_read));
  Fq_read  *seq2 = malloc(sizeof(Fq_read));
  DS_read *ptr_dsread = malloc(sizeof(DS_read));
  init_dsLUTs();
  char r1a[60], r2a[60];
  do {
    newl1 = fread(buffer1+offset1, 1, B_LEN-offset1, f_r1);
    newl2 = fread(buffer2+offset2, 1, B_LEN-offset2, f_r2);
     newl1 += offset1; 
     newl2 += offset2; 
     buffer1[newl1] = '\0';
     buffer2[newl2] = '\0';
     while (buffer1[j1] != 0 ||  buffer2[j2] != 0) {
        if (buffer1[j1] == '\n') {
          l1_f = j1;
          get_fqread(seq1, buffer1, l1_i, l1_f, nl1, par_TF.L, 0);
          if ((nl1++ % 4 == 3)) {
             stop1 = 1;
             j1++;
          }
          l1_i = l1_f + 1;
        }
        if (buffer2[j2] == '\n') {
          l2_f = j2; 
          get_fqread(seq2, buffer2, l2_i, l2_f, nl2, par_TF.L, 0); 
          if ((nl2++ % 4 == 3)) {
             stop2 = 1;
             j2++;
          }
          l2_i = l2_f + 1; 
        }
        if (stop1 && stop2 ) { // Do the stufff!! 
           nentries++;
           match_dsread(seq1, seq2, ptr_dsread);
           strncpy(r1a,(char *)seq1 -> line2, ptr_dsread -> posmax - 1);
           strncpy(r2a,(char *)seq2 -> line2, ptr_dsread -> posmax - 1);
           r1a[ptr_dsread-> posmax-1] = '\0'; 
           r2a[ptr_dsread-> posmax-1] = '\0';
           printf( "\t\t\t\t\t\t\t READ N %d  \n", nentries);
           printf("READ1: %s || ", seq1 -> line2);
           printf("READ2: %s\n", seq2 -> line2);
           if(ptr_dsread -> nmax > 60) {
             printf("READ1: %s|", r1a);
             printf("\033[32;1m%s\033[0m || ", seq1 -> line2 + ptr_dsread -> posmax);
             printf("READ2: %s|", r2a);
             printf("\033[32;1m%s\033[0m ", seq2 -> line2 + ptr_dsread -> posmax );
             printf("Match in pos: %d\n", ptr_dsread -> posmax);
           }  else {
           printf( "\t\t\t\t\t\t\t NO MATCH FOUND  \n");
           }
           // reset stop
           stop1 = 0; 
           stop2 = 0; 
        } else if (stop1 && !stop2 ) {
           j2++;
        } else if (!stop1 && stop2) {
           j1++;
        } else {
           j1++; 
           j2++;   
        }
     }  // end buffer loop
     offset1 = newl1 - l1_i;
     if (offset1 > -1)
       memcpy(buffer1, buffer1 + l1_i, offset1);
     l1_f = -1;
     l1_i = 0;
     offset2 = newl1 - l2_i;
     if (offset2 > -1)
       memcpy(buffer2, buffer2 + l2_i, offset2);
     l2_f = -1;
     l2_i = 0;
  }  while ((newl1 > offset1) || (newl2 > offset2));  // end read buffer
  fclose(f_r1);
  fclose(f_r2);
  return 0; 
}
