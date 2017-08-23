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
 * @file init_Sreport.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 09.08.2017
 * @brief Help dialog for Sreport and initialization of
 * the command line arguments.
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "init_Sreport.h"
#include "config.h"

extern Iparam_Sreport par_SR;

/**
 * @brief Function that prints Sreport help dialog when called.
*/
void printHelpDialog_Sreport() {
  const char dialog[] =
    "Usage: ./Sreport -i <INPUT_FOLDER> -o <OUTPUT_FILE> \n"
    "Uses all *bin files found in a folder ( output of Qreport) \n"
    "and generates a summary report in html format.\n"
    "Options:\n"
     " -v Prints package version.\n"
     " -h Prints help dialog.\n"
     " -i Input folder containing *bin data (output from Qreport)."
     " Mandatory option.\n"
     " -o Output file (with NO extension). Mandatory option.\n\n";
  fprintf(stderr, "%s", dialog);
}

/**
 * @brief Reads in the arguments passed through the command line to Sreport.
 *   and stores them in the global variable par_SR.
 *
*/
void getarg_Sreport(int argc, char **argv) {
  if (argc != 2 && argc != 5) {
     fprintf(stderr, "Not adequate number of arguments\n");
     printHelpDialog_Sreport();
     fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
  }
  char option;
  while ((option = getopt(argc, argv, "hvi:l:t:q:n:o:f:")) != -1) {
    switch (option) {
      case 'h':  // show the HelpDialog
        printHelpDialog_Sreport();
        exit(EXIT_SUCCESS);
        break;
      case 'v':  // Print version
        printf("Sreport version %s \nWritten by Paula Perez Rubio\n", VERSION);
        exit(EXIT_SUCCESS);
        break;
      case 'i':
        par_SR.inputfolder = optarg;
      case 'o':
        snprintf(par_SR.outputfile, MAX_FILENAME, "%s.html" , optarg);
      default:
        fprintf(stderr, "%s: option `-%c' is invalid: ignored\n",
                              argv[0], option);
        printHelpDialog_Sreport();
        fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
        exit(EXIT_FAILURE);
        break;
    }
  }
  if (par_SR.inputfolder == NULL) {
     fprintf(stderr, "Input folder was not properly initialized. \n");
     fprintf(stderr, "Exiting program.\n");
     printHelpDialog_Sreport();
     fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
  }
  if (par_SR.outputfile == NULL) {
     fprintf(stderr, "html output file was not properly initialized. \n");
     fprintf(stderr, "Exiting program.\n");
     printHelpDialog_Sreport();
     fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
     exit(EXIT_FAILURE);
  }
}

