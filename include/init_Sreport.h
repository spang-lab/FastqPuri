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
 * @file init_Sreport.h 
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 09.08.2017
 * @brief Help dialog for Sreport and initialization of 
 * the command line arguments.
 */


#ifndef INIT_SREPORT_H_
#define INIT_SREPORT_H_

#include "defines.h"

/**
 * @brief contains Sreport input parameters 
 */
typedef struct _iparam_Sreport{
  char *inputfolder; /**< input folder */
  char outputfile[MAX_FILENAME]; /**< html outputfile path */
  char *Rmd_file; /**< Rmd file path */ 
  char pBuf[MAX_FILENAME]; /**< html outputfile path */
} Iparam_Sreport; 

void printHelpDialog_Sreport();

void getarg_Sreport(int argc, char **argv);

#endif  // endif INIT_SREPORT_H_
