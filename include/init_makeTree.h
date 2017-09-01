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
 * @file init_makeTree.h
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 23.08.2017
 * @brief Help dialog for makeTree and initialization of
 * the command line arguments.
 */

#ifndef INIT_MAKETREE_H_
#define INIT_MAKETREE_H_

#include "defines.h"

/**
 * @brief contains makeTree input parameters
 */
typedef struct _iparam_makeTree {
  char *inputfasta; /**< fasta input file */
  char outputfile[MAX_FILENAME]; /**< outputfile path */
  int L; /**< tree depth */
} Iparam_makeTree;

void printHelpDialog_makeTree();

void getarg_makeTree(int argc, char **argv);

#endif  // endif INIT_MAKETREE_H_
