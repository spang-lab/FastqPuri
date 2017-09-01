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
 * @file Rcommand_Sreport.c
 * @brief get Rscript command for Sreport
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 09.08.2017
 *
 */

#include <stdio.h>
#include <unistd.h>
#include "Rcommand_Sreport.h"
#include "init_Sreport.h"
#include "defines.h"
#include "config.h"

extern Iparam_Sreport par_SR; /**< input parameters Sreport*/

/**
 * @brief returns Rscript command that generates the summary report in html
 * @code{.R}
 * # To run between quotation marks after:  Rscript_RBioC -e (Rscript)
 * inputfolder = normalizePath( <par.SR.inputfolder>, mustWork = TRUE);
 * output = <par_SR.outputfile>;
 * output_file = gsub('.* /', '', output);
 * path = gsub('[^/]+$', '', output);
 * if (path != '') {
 *   outputfile = paste0(normalizePath(path, mustWork = TRUE), '/', outputfile);
 * } else {
 *   outputfile = paste0(cwd, '/', output_file); # cwd: current working dir
 * }; 
 * rmarkdown::render(<par_SR.Rmd_file>, 
 *                   params = list(inputfolder = inputfolder, version= VERSION),
 *                   output_file = output_file)
 * @endcode
 * */
char *command_Sreport(){
  char command[MAX_RCOMMAND];
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL)
      fprintf(stdout, "Current working dir: %s\n", cwd);
  else
      perror("getcwd() error");
  snprintf(command, MAX_RCOMMAND, "%s -e \" inputfolder = normalizePath('%s', \
mustWork = TRUE); output = '%s';\
output_file = gsub('.*/', '', output);\
path = gsub('[^/]+$', '', output);\
if (path != '') { output_file = paste0(normalizePath(path, mustWork=TRUE)\
, '/', output_file); \
} else {\
output_file = paste0('%s', '/', output_file); };\
rmarkdown::render('%s', params = list(inputfolder = inputfolder, \
version = '%s'), output_file = output_file)\"",
        RSCRIPT_EXEC, par_SR.inputfolder,
        par_SR.outputfile, cwd, par_SR.Rmd_file, VERSION);
  char *str_command = command;
  return str_command;
}
