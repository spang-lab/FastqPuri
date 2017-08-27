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
 * @file init_treefilter.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 24.08.2017
 * @brief help dialog for trimFilter and initialization of the 
 * command line arguments.
 *
 * */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include "init_trimFilter.h"
#include "defines.h"
#include "config.h"

extern Iparam_trimFilter par_TF; /**< Input parameters of makeTree */

/**
 * @brief Function that prints trimFilter help dialog when called.
*/
void printHelpDialog_trimFilter() {
  const char dialog[] = 
   "Usage: trimFilter --ifq <INPUT_FILE.fq> --length <READ_LENGTH> \n"
   "                  --output [O_PREFIX] --minQ [MINQ]  --idx [INDEX_FILE] \n"
   "                  --ifa [INPUT.fa] --method [TREE|SA|BLOOM] \n"
   "                  --trimQ [NO|ALL|ENDS|FRAC|ENDSFRAC|GLOBAL] | \n"
   "                  (--percent percent) | (--global n1 n2) | (--minL minL)\n"
   "                  --trimN [NO|ALL|ENDS|STRIP]  \n"
   "Reads in a fq file (gz, bz2, z formats also accepted) and removes: \n"
   "  * low quality reads,\n"
   "  * reads containing N base callings,\n"
   "  * reads representing contaminations, belonging to sequences in INPUT.fa\n"
   "Outputs 4 [O_PREFIX]_fq.gz files containing: \"good\" reads, discarded \n"
   "low Q reads discarded reads containing N's, discarded contaminations.\n"
   "Options:\n"
   " -v, --version prints package version.\n"
   " -h, --help    prints help dialog.\n"
   " -f, --ifq     fastq input file [*fq|*fq.gz|*fq.bz2], mandatory option.\n"
   " -l, --length  read length: length of the reads, mandatory option.\n"
   " -o, --output  output prefix (with path), optional (default ./out).\n"
   " -q, --minQ    minimum quality allowed (int), optional (default 27).\n"
   " -x, --idx     index input file: optional.\n"
   " -a, --ifa     fasta input file [*fa|*fa.gz|*fa.bz2]. It only makes sense\n"
   "               when the method TREE is used. Not allowed otherwise.\n"
   " -k, --kmer    kmer length: Length of the Kmers to be searched in the \n"
   "               contaminations file *fa. If an index file is given, the \n"
   "               length specified in it is used.\n"
   " -M, --method  method used to look for contaminations: \n"
   "               TREE:  uses a 4-ary tree. Index file optional,\n"
   "               SA:    uses a suffix array. Index file mandatory,\n"
   "               BLOOM: uses a bloom filter. Index file mandatory.\n"
   " -Q, --trimQ   NO:       does nothing to low quality reads (default),\n"
   "               ALL:      removes all reads containing at least one low\n"
   "                         quality nucleotide.\n" 
   "               ENDS:     trims the ends of the read if their quality is\n"
   "                         below the threshold -q,\n"
   "               FRAC:     discards a read if the fraction of bases whose\n"
   "                         quality lies below \n" 
   "                         the threshold -q is over 5\% or a user defined\n"
   "                         percentage in -p.\n"
   "               ENDSFRAC: trims the ends and then discards the read if \n" 
   "                         there are more low quality nucleotides than the\n"
   "                         allowed by the option -p.\n"
   "               GLOBAL:   removes n1 cycles on the left and n2 on the \n"
   "                         right, specified in -g.\n"
   "               All reads are discarded if they are shorter than minL.\n"
   " -m, --minL    minimum length allowed for a read before it is discarded\n"
   "               (default 25).\n"
   " -p, --percent percentage of low quality bases to be admitted before \n"
   "               discarding a read (default 5), \n" 
   " -g, --global  required option if --trimQ GLOBAL is passed. Two int,\n"
   "               n1, n2, have to be passed specifying the number of cycles \n"
   "               to be globally cut from the left and right, respectively.\n" 
   " -N, --trimN   NO:     does nothing to reads containing N's,\n" 
   "               ALL:    removes all reads containing N's,\n"
   "               ENDS:   trims ends of reads with N's,\n"
   "               STRIPS: looks for the largest substring with no N's.\n"
   "               All reads are discarded if they are shorter than minL.\n";
  fprintf(stderr, "%s", dialog);
}

/**
 * @brief Reads in the arguments passed through the command line to trimFilter.
 *   and stores them in the global variable par_TF.
*/
void getarg_trimFilter(int argc, char **argv) {
  if ( argc != 2 && (argc > 21 || argc % 2 == 0 || argc == 1) ) {
     fprintf(stderr, "Not adequate number of arguments");
     printHelpDialog_trimFilter();
     fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
  }
  static struct option long_options[] = { 
     {"version", no_argument, 0, 'v'},
     {"help", no_argument, 0, 'h'},
     {"ifq", required_argument, 0, 'f'},
     {"length", required_argument, 0, 'l'},
     {"output", required_argument, 0, 'o'},
     {"minQ", required_argument, 0, 'q'},
     {"idx", required_argument, 0, 'x'},
     {"ifa", required_argument, 0, 'a'},
     {"method", required_argument, 0, 'M'},
     {"trimQ", required_argument, 0, 'Q'},
     {"percent", required_argument, 0, 'm'},
     {"global", required_argument, 0, 'p'},
     {"minL", required_argument, 0, 'g'},
     {"trimN", required_argument, 0, 'N'},
  };
  int option; 
  int method_len = 20; 
  while ((option = getopt_long(argc, argv, "hvf:l:o:q:x:a:M:Q:m:p:g::N:", 
        long_options, 0)) != -1) {
    switch (option) {
      case 'h':
        printHelpDialog_trimFilter();
        exit(EXIT_SUCCESS);
        break; 
      case 'v':
        printf("Qreport version %s \nWritten by Paula Perez Rubio\n", VERSION);
        exit(EXIT_SUCCESS);
        break;
      case 'f':
         par_TF.Ifq = optarg; 
         break; 
      case 'l':
         par_TF.L = atoi(optarg);
         break; 
      case 'o':
         par_TF.Oprefix = optarg; 
         break; 
      case 'q':
         par_TF.minQ = atoi(optarg);
         break; 
      case 'x':
         par_TF.Iidx = optarg; 
         par_TF.is_idx = true; 
         break; 
      case 'a':
         par_TF.Ifa = optarg; 
         par_TF.is_fa = true;     
         break; 
      case 'M':
         par_TF.method = (!strncmp(optarg, "TREE", method_len)) ? TREE :
            (!strncmp(optarg, "SA", method_len)) ? SA :
            (!strncmp(optarg, "BLOOM", method_len)) ? BLOOM : ERROR; 
         break; 
      case 'Q':
         par_TF.trimQ = (!strncmp(optarg, "NO", method_len)) ? NO : 
            (!strncmp(optarg, "FRAC", method_len)) ? FRAC :
            (!strncmp(optarg, "ENDS", method_len)) ? ENDS :
            (!strncmp(optarg, "ENDSFRAC", method_len)) ? ENDSFRAC :
            (!strncmp(optarg, "GLOBAL", method_len)) ? GLOBAL : ERROR;
         break; 
      case 'm':
         par_TF.minL = atoi(optarg); 
         break; 
      case 'p':
         par_TF.percent = atoi(optarg);
         break; 
      case 'g':
         par_TF.globleft = atoi(optarg);
         par_TF.globright = atoi(argv[optind]);
         break; 
      case 'N':
         par_TF.trimN = (!strncmp(optarg, "NO", method_len)) ? NO : 
            (!strncmp(optarg, "ALL", method_len)) ? ALL :
            (!strncmp(optarg, "ENDS", method_len)) ? ENDS :
            (!strncmp(optarg, "STRIP", method_len)) ? STRIP : ERROR;
         break; 
      default: 
        fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                           argv[0], optopt);
        printHelpDialog_trimFilter();
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        break;
    }
  }
  fprintf(stderr, "Starting trim Filter.\n");
  // Checking the input
  // Ifq is a mandatory argument
  if (par_TF.Ifq == NULL) {
    printHelpDialog_trimFilter();
    fprintf(stderr, "Input *fq file name was not properly initialized and \n");
    fprintf(stderr, "is a mandatory option. (--ifq <INPUT_FILE.fq>)\n");
    fprintf(stderr, "Exiting program\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  } else {
    fprintf(stderr, "- Fastq input file: %s \n",par_TF.Ifq);
  }
  // Read length is a mandatory argument
  if (par_TF.L == 0) {
    printHelpDialog_trimFilter();
    fprintf(stderr, "Read length  was not properly initialized and \n");
    fprintf(stderr, "is a mandatory option. (--length <READ_LENGTH>)\n");
    fprintf(stderr, "Exiting program\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    exit(EXIT_FAILURE);
  } else {
    fprintf(stderr, "- Read length: %d.\n",par_TF.L);
  } 
  // handling output prefix
  if (par_TF.Oprefix == NULL) {
    par_TF.Oprefix = "./out"; 
    fprintf(stderr, "- Output prefix: %s (default value).\n", par_TF.Oprefix);
  } else {
    fprintf(stderr, "- Output prefix: %s.\n", par_TF.Oprefix);
  }
  // handling minQ 
  if (par_TF.minQ == 0) {
    par_TF.minQ = DEFAULT_MINQ; 
    fprintf(stderr, "- Min Quality: %d (default value).\n", par_TF.minQ);
  } else {
    fprintf(stderr, "- Min Quality: %d.\n", par_TF.minQ);
  }
  // Checking contamination search method 
  if (par_TF.method == 0) {
   fprintf(stderr, "- Not looking for contaminations\n");
   if (par_TF.is_fa || par_TF.is_idx) {
      fprintf(stderr, "OPTION_ERROR: No method specified but fasta or index\n");
      fprintf(stderr, "              file passed. Revise options (--help).\n");
      fprintf(stderr, "Exiting program\n");
      fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
  } else if (par_TF.method == TREE) {
     fprintf(stderr, "- Looking for contaminations with a tree struct.\n");
     if (par_TF.is_fa && par_TF.is_idx) {
        fprintf(stderr, "OPTION_ERROR: a fasta inputfile and an index input");
        fprintf(stderr, "file given\n");
        fprintf(stderr, "Mutually exclusive, revise options with --help\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     } else if (par_TF.is_fa) {
        fprintf(stderr, "- Constructing tree on the flight from %s.\n",
              par_TF.Ifa);
     } else if (par_TF.is_idx) {
        fprintf(stderr, "- Reading tree from file %s.\n",
              par_TF.Iidx);
     } else {
        fprintf(stderr, "OPTION_ERROR: a fasta file or an index file need \n");
        fprintf(stderr, "              to be specified. Revise options with ");
        fprintf(stderr, "--help.\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     }
  } else if (par_TF.method == SA) { 
     fprintf(stderr, "- Looking for contaminations with a suffix array.\n");
     if (par_TF.is_fa && par_TF.is_idx) {
        fprintf(stderr, "OPTION_ERROR: a fasta inputfile and an index input");
        fprintf(stderr, "file given\n");
        fprintf(stderr, "Mutually exclusive, revise options with --help\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     } else if (par_TF.is_fa) {
        fprintf(stderr, "OPTION_ERROR: fasta file passed as an argument. \n");
        fprintf(stderr, "              SA not computed on the flight.\n");
        fprintf(stderr, "              an index file needs to be provided.\n");
        fprintf(stderr, "              Revise your options with --help.\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     } else if (par_TF.is_idx) {
        fprintf(stderr, "- Reading SA from file %s.\n",
              par_TF.Iidx);
     } else {
        fprintf(stderr, "OPTION_ERROR: an index file needs ");
        fprintf(stderr, "to be specified. Revise options with ");
        fprintf(stderr, "--help.\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     }
  } else if (par_TF.method == BLOOM) { 
     fprintf(stderr, "- Looking for contaminations with a bloom filter.\n");
     if (par_TF.is_fa && par_TF.is_idx) {
        fprintf(stderr, "OPTION_ERROR: a fasta inputfile and an index input");
        fprintf(stderr, "file given\n");
        fprintf(stderr, "Mutually exclusive, revise options with --help\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     } else if (par_TF.is_fa) {
        fprintf(stderr, "OPTION_ERROR: fasta file passed as an argument. \n");
        fprintf(stderr, "              BF not computed on the flight.\n");
        fprintf(stderr, "              an index file needs to be provided.\n");
        fprintf(stderr, "              Revise your options with --help.\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     } else if (par_TF.is_idx) {
        fprintf(stderr, "- Reading Bloom filter from file %s.\n",
              par_TF.Iidx);
     } else {
        fprintf(stderr, "OPTION_ERROR: an index file needs ");
        fprintf(stderr, "to be specified. Revise options with ");
        fprintf(stderr, "--help.\n");
        fprintf(stderr, "Exiting program\n");
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
     }
  } else{
     fprintf(stderr, "OPTION_ERROR: Invalid --method option.\n");
     fprintf(stderr, "              Possible options: TREE, SA, BLOOM\n");
     fprintf(stderr, "              Revise your options with --help.\n");
     fprintf(stderr, "Exiting program\n");
     fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
  }
  // Check trimQ method 
  if (par_TF.trimQ == NO) {
    fprintf(stderr, "- Trimming low Q bases method: NO\n"); 
  } else if (par_TF.trimQ == ALL) {
    fprintf(stderr, "- Trimming low Q bases method: ALL\n"); 
  } else if (par_TF.trimQ == ENDS) {
    fprintf(stderr, "- Trimming low Q bases method: ENDS\n"); 
  } else if (par_TF.trimQ == FRAC) {
    fprintf(stderr, "- Trimming low Q bases method: FRAC\n"); 
    if (par_TF.percent == 0) {
      par_TF.percent = 5; 
    } 
    par_TF.nlowQ = par_TF.L*par_TF.percent/100+1; 
    fprintf(stderr, "- Read discarded if containing more than %d %c",
           par_TF.percent,'%');
    fprintf(stderr,"lowQ bases (< %d).\n", par_TF.nlowQ); 
  } else if (par_TF.trimQ == ENDSFRAC) {
    if (par_TF.percent == 0) {
      par_TF.percent = 5; 
    } 
    par_TF.nlowQ = par_TF.L*par_TF.percent/100+1; 
    fprintf(stderr, "- Trimming low Q bases method: ENDSFRAC\n");
    fprintf(stderr, "- Trimmed read discarded if containing more than %d %c",
           par_TF.percent,'%');
    fprintf(stderr,"lowQ bases (< %d).\n", par_TF.nlowQ); 
  } else if (par_TF.trimQ == GLOBAL) {
    fprintf(stderr, "- Trimming low Q bases method: GLOBAL\n");
    fprintf(stderr, "- Trimming globally %d from left and %d from right\n", 
           par_TF.globleft, par_TF.globright);
  } else {
      fprintf(stderr, "OPTION_ERROR: Invalid --trimQ option.\n");
      fprintf(stderr, "              Possible options: NO, ALL, ENDS,");
      fprintf(stderr, " ENDSFRAC, GLOBAL.\n");
      fprintf(stderr, "              Revise your options with --help.\n");
      fprintf(stderr, "Exiting program\n");
      fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
  }
  // Consistenty checks 
  if ((par_TF.trimQ != ENDS) && (par_TF.trimQ != ENDSFRAC) && 
      (par_TF.percent != 0)) { 
      fprintf(stderr, "OPTION_ERROR: --percent passed as an option, but ");
      fprintf(stderr, "neither ENDS nor\n");
      fprintf(stderr, "              ENDSFRAC where passed to --trimQ. \n");
      fprintf(stderr, "              Maybe you meant something else?. \n");
      fprintf(stderr, "              Revise your options with --help.\n");
      fprintf(stderr, "Exiting program\n");
      fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
  if ((par_TF.trimQ != GLOBAL) && ((par_TF.globleft != 0) || 
           (par_TF.globright != 0) )) { 
      fprintf(stderr, "OPTION_ERROR: --global passed as an option, but ");
      fprintf(stderr, "GLOBAL not passed to --trimQ\n");
      fprintf(stderr, "              Maybe you meant something else?. \n");
      fprintf(stderr, "              Revise your options with --help.\n");
      fprintf(stderr, "Exiting program\n");
      fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }
   if (par_TF.trimN == NO) {
      fprintf(stderr, "- Trimming reads with N's, method: NO\n"); 
   } else if (par_TF.trimN == ALL) {
      fprintf(stderr, "- Trimming reads with N's, method: ALL\n"); 
   } else if (par_TF.trimN == ENDS) {
      fprintf(stderr, "- Trimming reads with N's, method: ENDS\n"); 
   } else if (par_TF.trimN == STRIP) {
      fprintf(stderr, "- Trimming reads with N's, method: STRIP\n"); 
   } else {
      fprintf(stderr, "OPTION_ERROR: Invalid --trimN option.\n");
      fprintf(stderr, "              Possible options: NO, ALL, ENDS,");
      fprintf(stderr, " STRIP.\n");
      fprintf(stderr, "              Revise your options with --help.\n");
      fprintf(stderr, "Exiting program\n");
      fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
      exit(EXIT_FAILURE);
   }

   if (par_TF.minL == 0) {
     par_TF.minL = DEFAULT_MINL;
     fprintf(stderr, "- Minimum accepted trimmed read length: %d (default).\n",
            par_TF.minL); 
   } else {
     fprintf(stderr, "- Minimum accepted trimmed read length: %d.\n",
            par_TF.minL); 
   }
}
