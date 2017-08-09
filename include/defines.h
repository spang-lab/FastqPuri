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
 * @file defines.h
 * @brief Macro definitions
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 07.08.2017
 * 
 */

#ifndef DEFINES_H_
#define DEFINES_H_


#define MAX_FILENAME 300 /**< Maximum # chars in a filename*/
#define DEFAULT_MINQ 27   /**< Minimum quality threshold */ 
#define DEFAULT_NTILES 96 /**< Default number of tiles */ 
#define DEFAULT_NQ 46 /**< Default number of different quality values */ 
#define ZEROQ 33 /**< ASCII code of lowest quality value (!) */ 
#define N_ACGT 5 /**< Number of different nucleotides in the fq file */ 
#define MAX_RCOMMAND  4000 /**< Maximum # chars in R command*/


#endif

