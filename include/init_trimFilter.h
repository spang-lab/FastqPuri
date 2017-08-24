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
 * @file init_treefilter.h
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 24.08.2017
 * @brief help dialog for trimFilter and initialization of the 
 * command line arguments.
 *
 * */

#ifndef INIT_TRIMFILTER_H
#define INIT_TRIMFILTER_H

#include <stdio.h>
#include "defines.h" 

typedef struct _iparam_trimFilter {
  char *Ifq;
  char *Ifa;
  char *Iidx;
  char *Oprefix;
  int trimQ;  // NO(0), FRAC(1), ENDS(2), ENDSFRAC(3), GLOBAL(4)
  int trimN;  // NO(0), ALL(1), ENDS(2), STRIP(3)
  bool is_fa, is_idx, tree;
  double score; 
  int  minQ, L, minL, nlowQ, Lmer_len, globleft, globright, percent;
} Iparam_trimFilter;

void printHelpDialog_trimFilter();

void getarg_Qreport(int argc, char **argv);

#endif  // INIT_TRIMFILTER_H
