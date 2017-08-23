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

// General 
#define B_LEN 131072  /**< buffer size */
#define MAX_FILENAME 300 /**< Maximum # chars in a filename */
#define bool short  /**< define a bool type */
#define true 1   /**< assign true to 1 */
#define false 0  /**< assign false to 0 */

#ifndef max
   #define max( a, b ) ( ((a) > (b)) ? (a) : (b) ) /**< max function */
#endif

#ifndef min
   #define min( a, b ) ( ((a) < (b)) ? (a) : (b) ) /**< min function */
#endif

#ifndef mem_usageMB
   #define mem_usageMB()  fprintf(stderr, \
         "- Current allocated memory: %ld MB.\n", \
         alloc_mem >> 20)  /**< returns allocated memory in MB */
#endif

#ifndef mem_usage
   #define mem_usage()  fprintf(stderr, \
         "- Current allocated memory: %ld Bytes.\n", \
         alloc_mem) /**< returns allocated memory in Bytes */
#endif


// Q_report, S_report
#define DEFAULT_MINQ 27   /**< Minimum quality threshold */ 
#define DEFAULT_NTILES 96 /**< Default number of tiles */ 
#define DEFAULT_NQ 46 /**< Default number of different quality values */ 
#define ZEROQ 33 /**< ASCII code of lowest quality value (!) */ 
#define N_ACGT 5 /**< Number of different nucleotides in the fq file */ 
#define MAX_RCOMMAND  4000 /**< Maximum # chars in R command*/


// Fasta files
#define FA_ENTRY_BUF 20 /**< buffer for fasta entries*/

// Tree 
#define T_ACGT 4 /**< Number of children per node in tree*/
#define NPOOL_1D 1048576 /**< Number of Node structs allocated in inner dim */
#define NPOOL_2D 16  /**< Number of *Node allocated in outer dim */
#define MAX_FASZ_TREE 1e7 /**< Maximum fasta size for constructing a tree.
                               DECIDE A SENSIBLE SIZE! */



#endif

