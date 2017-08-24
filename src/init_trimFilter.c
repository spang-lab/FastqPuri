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

#include "init_trimFilter.h"

extern Iparam_trimFilter par_TF; /**< Input parameters of makeTree */

/**
 * @brief Function that prints trimFilter help dialog when called.
*/
void printHelpDialog_trimFilter() {
  const char dialog[] = 
   "Usage ./trimFilter ";

  fprintf(stderr, "%s", dialog);
}

/**
 * @brief Reads in the arguments passed through the command line to trimFilter.
 *   and stores them in the global variable par_TF.
*/
void getarg_trimFilter(int argc, char **argv) {

}
