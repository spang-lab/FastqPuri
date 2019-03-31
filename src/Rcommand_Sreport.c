/****************************************************************************
 * Copyright (C) 2017 by Paula Perez Rubio                                  *
 *                                                                          *
 * This file is part of FastqPuri.                                      *
 *                                                                          *
 *   FastqPuri is free software: you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as                *
 *   published by the Free Software Foundation, either version 3 of the     *
 *   License, or (at your option) any later version.                        *
 *                                                                          *
 *   FastqPuri is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU General Public License for more details.                           *
 *                                                                          *
 *   You should have received a copy of the GNU General Public License      *
 *   along with FastqPuri.                                              *
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
#include <stdlib.h>
#include <libgen.h>
#include <strings.h>
#include <unistd.h>
#include "tinydir.h"
#include "Rcommand_Sreport.h"
#include "init_Sreport.h"
#include "copy_file.h"
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
char *command_Sreport(char **new_dir_ptr){
  char *command = calloc(MAX_RCOMMAND, sizeof(char));
  command[0] = '\0';
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL)
      fprintf(stderr, "- Current working dir: %s\n", cwd);
  else
      perror("getcwd() error");
  tinydir_file file;
  tinydir_dir dir;
  tinydir_open(&dir, par_SR.inputfolder);
  int hasBin = 0;
  while (dir.has_next) {
    tinydir_readfile(&dir, &file);
    if (strcmp(file.extension, "bin") == 0) {
      hasBin = 1;
      break;
    }
    tinydir_next(&dir);
  }
  tinydir_close(&dir);
  if (hasBin == 0) {
    fprintf(stderr, "Sreport did not find any binary files in %s.\n", par_SR.inputfolder);
    fprintf(stderr, "Maybe something went wrong with Qreport trying to generate such files.\n");
  } else {
#ifdef HAVE_RPKG
  char template[] = "/tmp/FastqPuri_XXXXXX";
  char *new_dir = mkdtemp(template);
  *new_dir_ptr = new_dir;
  char old_dir_tmp[MAX_FILENAME];
  strncpy(old_dir_tmp, par_SR.Rmd_file, MAX_FILENAME-1);
  char *old_dir = dirname(old_dir_tmp);
  char rmd_summary_report_name_tmp[MAX_FILENAME];
  strncpy(rmd_summary_report_name_tmp, par_SR.Rmd_file, MAX_FILENAME-1);
  char *rmd_summary_report_name = basename(rmd_summary_report_name_tmp);
  char style_fname_old[MAX_FILENAME], utils_fname_old[MAX_FILENAME];
  char style_fname_new[MAX_FILENAME], utils_fname_new[MAX_FILENAME];
  char rmd_summary_report_new[MAX_FILENAME];
 

  snprintf(rmd_summary_report_new, MAX_FILENAME, "%s/%s", new_dir, rmd_summary_report_name);  
  snprintf(style_fname_old, MAX_FILENAME, "%s/style.css", old_dir);  
  snprintf(style_fname_new, MAX_FILENAME, "%s/style.css", new_dir);  
  snprintf(utils_fname_old, MAX_FILENAME, "%s/utils.R", old_dir);  
  snprintf(utils_fname_new, MAX_FILENAME, "%s/utils.R", new_dir);  
   
  copy_file(par_SR.Rmd_file, rmd_summary_report_new);
  copy_file(utils_fname_old, utils_fname_new);
  copy_file(style_fname_old, style_fname_new);
 
   
  
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
        par_SR.outputfile, cwd, rmd_summary_report_new, VERSION);
#endif
  }
  return command;
}
