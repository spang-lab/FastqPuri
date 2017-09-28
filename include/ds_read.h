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
 * @file ds_read.h
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 20.09.2017
 * @brief fastq entries manipulations (read/write)
 *
 * */

#include "fq_read.h"
#include "config.h"

#define COMPACT 2
#define READ_HALFLEN READ_MAXLEN/2

/*
 * @brief  stores two paired fastq entries, r1, r2.
 *
 * */
typedef struct ds_read {
   Fq_read *r1;  /** Fq_read struct containing read1  */
   Fq_read *r2;  /** Fq_read struct containing read2  */
   int Lhalf1;    /** length of the compactified strings read1*/ 
   int Lhalf2;    /** length of the compactified strings read2*/ 
   unsigned char c1[READ_HALFLEN]; /** compact string encoding read1 */
   unsigned char c1sh[READ_HALFLEN]; /** compact string encoding read1 shifted*/
   unsigned char c2[READ_HALFLEN]; /** compact string encoding read2 */ 
   int nmax; /** Maximum number of matches*/ 
   int posmax; /** Position of palindromic match*/
} DS_read;

void init_dsLUTs(); 
void match_dsread(Fq_read *r1, Fq_read *r2, DS_read *ptr_dsread);
