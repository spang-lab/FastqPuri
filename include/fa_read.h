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
 * @file fa_read.h
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 16.08.2017
 * @brief reads in and stores fasta files
 *
 * */


#ifndef FA_READ_H
#define FA_READ_H

#include <stdint.h>

/**
 * @brief fasta entry
 * */
typedef struct _fa_entry {
  uint64_t  N;        /**<  Entry length (chars)*/
  char *seq; /**< sequence */
} Fa_entry;

/**
 * @brief stores sequences of a fasta file
 * */
typedef struct _fa_data {
  uint64_t nlines;    /**< Number of lines  in *fa file */
  int nentries;       /**< Number of entries in *fa file */
  int linelen;        /**< Line length of the *fa file entries */
  uint64_t *entrylen; /**< Array containing the length of the entries */
  Fa_entry *entry;  /**< Array with fasta entries (see Fa_entry)*/
} Fa_data;


int read_fasta(char *filename, Fa_data *ptr_fa);
void free_fasta(Fa_data *ptr_fa);

#endif
