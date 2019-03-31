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
 * @file Rcommand_Qreport.c
 * @brief get Rscript command for Qreport
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 07.08.2017
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <libgen.h>
#include <unistd.h>
#include <string.h>
#include "Rcommand_Qreport.h"
#include "init_Qreport.h"
#include "copy_file.h"
#include "defines.h"
#include "config.h"

extern Iparam_Qreport par_QR; /**< input parameters Qreport*/

/**
 * @brief returns Rscript command that generates the quality report in html
 * */
char *command_Qreport(char ** new_dir_ptr) {
  char *command = calloc(MAX_RCOMMAND,sizeof(char));
  char cwd[1024];
  if (getcwd(cwd, sizeof(cwd)) != NULL)
      fprintf(stderr, "- Current working dir: %s\n", cwd);
  else
      perror("getcwd() error");
#ifdef HAVE_RPKG
  // find the calling program to deduce where to find Rmd-files
  char szTmp[32]; char pBuf[MAX_FILENAME];
  size_t len = sizeof(pBuf);
  sprintf(szTmp, "/proc/%d/exe", getpid());
  int bytes = readlink(szTmp, pBuf, len);
  if ((size_t)bytes > len-1) bytes = len-1;
  if(bytes == 0) {
    fprintf(stderr, "Unexpected error when searching for call directoy!\n");
    exit(1);
  }
  char *old_dir = dirname(pBuf);
  if (strcmp(old_dir, INSTALL_DIR) != 0) {
    old_dir = dirname(old_dir);
    strcat(old_dir, "/R");
  } else {
    strcpy(pBuf, RMD_QUALITY_REPORT);
    old_dir = dirname(pBuf);
  }
//  fprintf(stderr, ">>>old_dir: %s, INSTALL_DIR: %s,\n   RSCRIPT_EXEC: %s, CMAKE_INSTALL_PREFIX: %s\n", old_dir, INSTALL_DIR, RSCRIPT_EXEC, CMAKE_INSTALL_PREFIX);
  
  char template[] = "/tmp/FastqPuri_XXXXXX";
//  fprintf(stderr, ">>>template: %s\n", template);
  char *new_dir = mkdtemp(template);
//  fprintf(stderr, ">>>template after mkdtemp: %s\n", template);
  *new_dir_ptr = new_dir;
//  fprintf(stderr, ">>>new_dir: %s\n", new_dir);
    
  char rmd_quality_report_name_tmp[] = RMD_QUALITY_REPORT;
  char *rmd_quality_report_name = basename(rmd_quality_report_name_tmp);
  
  char style_fname_old[MAX_FILENAME], utils_fname_old[MAX_FILENAME];
  char style_fname_new[MAX_FILENAME], utils_fname_new[MAX_FILENAME];
  char rmd_quality_report_old[MAX_FILENAME]; char rmd_quality_report_new[MAX_FILENAME];
  
  snprintf(rmd_quality_report_old, MAX_FILENAME, "%s/%s", old_dir, rmd_quality_report_name);  
  snprintf(rmd_quality_report_new, MAX_FILENAME, "%s/%s", new_dir, rmd_quality_report_name);  
//  fprintf(stderr, ">>>rmd_quality_report: %s -> %s\n", rmd_quality_report_old, rmd_quality_report_new);
  snprintf(style_fname_old, MAX_FILENAME, "%s/style.css", old_dir);  
  snprintf(style_fname_new, MAX_FILENAME, "%s/style.css", new_dir);  
//  fprintf(stderr, ">>>style_fname: %s -> %s\n", style_fname_old, style_fname_new);
  snprintf(utils_fname_old, MAX_FILENAME, "%s/utils.R", old_dir);  
  snprintf(utils_fname_new, MAX_FILENAME, "%s/utils.R", new_dir);  
//  fprintf(stderr, ">>>utils_fname: %s\n", utils_fname_old, utils_fname_new);
  fprintf(stderr, "- Rmd-file used to generate HTML: %s\n", rmd_quality_report_old);

  copy_file(rmd_quality_report_old, rmd_quality_report_new);
  copy_file(utils_fname_old, utils_fname_new);
  copy_file(style_fname_old, style_fname_new);

  snprintf(command, MAX_RCOMMAND,"%s -e \" inputfile = normalizePath('%s', \
 mustWork = TRUE) ; output = '%s';\
output_file = gsub('.*/', '', output);\
path = gsub('[^/]+$', '', output);\
if (path != '') { output_file = paste0(normalizePath(path, mustWork=TRUE)\
, '/', output_file); \
} else {\
output_file = paste0('%s', '/', output_file); };\
rmarkdown::render('%s', params = list(inputfile = inputfile, filter=%d, \
 version = '%s'), output_file = output_file)\"", RSCRIPT_EXEC,
       par_QR.outputfilebin, par_QR.outputfilehtml, cwd,
       rmd_quality_report_new, par_QR.filter, VERSION);
#endif
  return command;
}

