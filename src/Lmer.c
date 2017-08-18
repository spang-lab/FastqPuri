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
 * @file Lmer.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 18.08.2017
 * @brief Manipulation of Lmers and sequences
 *
 * */

#include "Lmer.h"
#include <stdio.h>

char LT[256]; /**< global variable. Lookup table. */

/**
 * @brief Initialize lookup table LT.
 *
 * {'a','c','g','t'}  --> {'\000','\001','\002','\003'}, rest '\004'.
 *
 * */
void init_map() {
  int i;
  for (i = 0; i < 256; i++)
     LT[i] = 4;
  LT['a'] = 0;
  LT['c'] = 1;
  LT['g'] = 2;
  LT['t'] = 3;
  LT['A'] = 0;
  LT['C'] = 1;
  LT['G'] = 2;
  LT['T'] = 3;
}

/**
 * @brief Initialize lookup table LT (for SA)
 *
 * {'a','c','g','t'}  --> {'\001','\002','\003','\004'}, rest '\005'.
 *
 * */
void init_map_SA() {
  int i;
  for (i = 0; i < 256; i++)
     LT[i] = 5;
  LT['a'] = 1;
  LT['c'] = 2;
  LT['g'] = 3;
  LT['t'] = 4;
  LT['A'] = 1;
  LT['C'] = 2;
  LT['G'] = 3;
  LT['T'] = 4;
}

/**
 * @brief Transforms an Lmer to the convention stored in the lookup table LT.
 *
 * */
void Lmer_sLmer(char* Lmer, int L) {
  int i;
  for (i = 0; i < L; i++){
     Lmer[i] = LT[(int)Lmer[i]];
  }
}

/**
 * @brief Obtains the reverse complement, for {'\000','\001','\002','\003'}.
 * */
void rev_comp(char *sLmer, int L) {
  char RC[5]= {3, 2, 1, 0, 4};
  int c, i, j;
  for (i = 0, j = L-1; i < j; i++, j--) {
     c = RC[(int) sLmer[i]];
     sLmer[i] = RC[(int) sLmer[j]];
     sLmer[j] = c;
  }
  if (i == j) {
     sLmer[i] = RC[(int) sLmer[j]];
  }
}

/**
 * @brief Obtains the reverse complement, for {'\001','\002','\003','\004'}.
 * */
void rev_comp2(char *sLmer, int L) {
  char RC[5]= {4, 3, 2, 1, 5};
  int c, i, j;
  for (i = 0, j = L-1; i < j; i++, j--) {
     c = RC[(int) sLmer[i] -1];
     sLmer[i] = RC[(int) sLmer[j]-1];
     sLmer[j] = c;
  }
  if (i == j) {
     sLmer[i] = RC[(int) sLmer[j]-1];
  }
}
