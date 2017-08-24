/****************************************************************************
 * copyright (c) 2017 by paula perez rubio                                  *
 *                                                                          *
 * this file is part of fastqarazketa.                                      *
 *                                                                          *
 *   fastqarazketa is free software: you can redistribute it and/or modify  *
 *   it under the terms of the gnu general public license as                *
 *   published by the free software foundation, either version 3 of the     *
 *   license, or (at your option) any later version.                        *
 *                                                                          *
 *   fastqarazketa is distributed in the hope that it will be useful,       *
 *   but without any warranty; without even the implied warranty of         *
 *   merchantability or fitness for a particular purpose.  see the          *
 *   gnu general public license for more details.                           *
 *                                                                          *
 *   you should have received a copy of the gnu general public license      *
 *   along with fastqarazketa.                                              *
 *   if not, see <http://www.gnu.org/licenses/>.                            *
 ****************************************************************************/

/**
 * @file trimFilter.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 25.08.2017
 * @brief trimFilter main function
 *
 * This file contains the trimFilter main function. BLA BLA 
 * See README_trimFilter.md for more details.
 *
 */

#include <stdio.h>
#include <time.h>
#include "init_trimFilter.h"
#include "init_Qreport.h"
#include "trim.h"


uint64_t alloc_mem = 0;  /**< global variable. Memory allocated in the heap.*/
Iparam_trimFilter par_TF;  /**< global variable: Input parameters of makeTree.*/
Iparam_Qreport par_QR;  /**< global variable: Input parameters of makeTree.*/

/**
 * @brief makeTree main function
 *
 * */
int main(int argc, char *argv[]) {
  clock_t start, end;
  double cpu_time_used;
  time_t rawtime;
  struct tm * timeinfo;

  // Start the clock
  start = clock();
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  
  // BODY of the function here! 


  // Obtaining elapsed time
  end = clock();
  cpu_time_used = (double)(end - start)/CLOCKS_PER_SEC;
  time(&rawtime);
  timeinfo = localtime(&rawtime);
  fprintf(stderr, "Finishing program at: %s", asctime(timeinfo) );
  fprintf(stderr, "Time elapsed: %f s.\n", cpu_time_used);
  return 0; 
}
