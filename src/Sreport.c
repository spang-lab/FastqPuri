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
 * @file Sreport.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 09.08.2017
 * @brief  Sreport main function
 *
 * This file contains the summary report main function. Given a folder
 * containing *bin as from Qreport output, Sreport generates a summary
 * report in html format. See README_Sreport.md for more details.
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "init_Sreport.h"
#include "Rcommand_Sreport.h"
#include "config.h"

Iparam_Sreport par_SR; /**< input parameters Sreport */


/**
 * @brief Qreport main function
 * */
int main(int argc, char *argv[]) {
  // Start the clock
  clock_t start, end;
  double cpu_time_used;
  time_t rawtime;
  struct tm * timeinfo;
  start = clock();
  time(&rawtime);
  timeinfo = localtime(&rawtime);

  // Get arguments
  getarg_Sreport(argc, argv);
  fprintf(stderr, "Starting program at: %s", asctime(timeinfo));
  fprintf(stderr, "- Input folder: %s\n", par_SR.inputfolder);
  fprintf(stderr, "- Output file: %s\n", par_SR.outputfile);
#ifdef HAVE_RPKG
  char * command = command_Sreport();
  printf("Running command: %s \n", command);
  int status;
  if ((status = system(command)) != 0) {
      fprintf(stderr, "Something went wrong when executing R script.\n");
      fprintf(stderr, "Most probably, a html file will not be generated.\n");
      fprintf(stderr, "Exiting program.\n");
      exit(EXIT_FAILURE);
  }
#else
  fprintf(stderr, "WARNING: html reports are NOT being generated.\n");
  fprintf(stderr, "         Dependencies not fulfilled.\n");
#endif

  // Obtaining elapsed time
  end = clock();
  cpu_time_used = (double)(end - start)/CLOCKS_PER_SEC;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  fprintf(stderr, "Finishing program at: %s", asctime(timeinfo) );
  fprintf(stderr, "Time elapsed: %f s.\n", cpu_time_used);

  return 0;
}
