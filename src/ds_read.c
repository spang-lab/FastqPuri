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
 * @file ds_read.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 20.09.2017
 * @brief fastq entries manipulations (read/write)
 *
 * */

#include <string.h>
#include <stdio.h>
#include "ds_read.h" 
#include "defines.h"

static uint8_t dsfw0[256], dsfw1[256], dsbw0[256], dsbw1[256];

void init_dsLUTs() {
  memset(dsfw0, 0x00, 256);
  memset(dsfw1, 0x00, 256);
  memset(dsbw0, 0x00, 256);
  memset(dsbw1, 0x00, 256);

 // Mappings
  dsfw0['a'] = 0x01; dsfw0['c'] = 0x02; dsfw0['g'] = 0x04; dsfw0['t'] = 0x08;
  dsfw0['A'] = 0x01; dsfw0['C'] = 0x02; dsfw0['G'] = 0x04; dsfw0['T'] = 0x08;
  dsfw1['a'] = 0x10; dsfw1['c'] = 0x20; dsfw1['g'] = 0x40; dsfw1['t'] = 0x80;
  dsfw1['A'] = 0x10; dsfw1['C'] = 0x20; dsfw1['G'] = 0x40; dsfw1['T'] = 0x80;

  dsbw0['a'] = 0x08; dsbw0['c'] = 0x04; dsbw0['g'] = 0x02; dsbw0['t'] = 0x01;
  dsbw0['A'] = 0x08; dsbw0['C'] = 0x04; dsbw0['G'] = 0x02; dsbw0['T'] = 0x01;
  dsbw1['a'] = 0x80; dsbw1['c'] = 0x40; dsbw1['g'] = 0x20; dsbw1['t'] = 0x10;
  dsbw1['A'] = 0x80; dsbw1['C'] = 0x40; dsbw1['G'] = 0x20; dsbw1['T'] = 0x10;
}

void match_dsread(Fq_read *r1, Fq_read *r2, DS_read *ptr_dsread) {
  ptr_dsread -> r1 = r1; 
  ptr_dsread -> r2 = r2; 
  // Length to be considered
  ptr_dsread -> Lhalf1 = (r1 -> L + 1)/2; 
  ptr_dsread -> Lhalf2 = (r2 -> L + 1)/2; 
  // Setting the data to 0
  memset(ptr_dsread -> c1, 0x00, READ_HALFLEN);
  memset(ptr_dsread -> c1sh, 0x00, READ_HALFLEN);
  memset(ptr_dsread -> c2, 0x00, READ_HALFLEN);
  ptr_dsread -> nmax = 0; 
  ptr_dsread -> posmax = 0; 
  // Creating the compact encoding
  int i, j = 0;

  for (i = 0; i < r1 -> L ; i++) {
    if(i%2) {
      ptr_dsread -> c1[j] |= dsfw1[(unsigned char)r1->line2[i]]; 
      if (i>0) 
        ptr_dsread -> c1sh[j] |= dsfw1[(unsigned char)r1->line2[i-1]]; 
      j++; 
    } else {
      ptr_dsread -> c1[j] |= dsfw0[(unsigned char)r1->line2[i]]; 
      ptr_dsread -> c1sh[j] |= dsfw0[(unsigned char)r1->line2[i-1]]; 
    }
  }
  j = 0; 
  for (i = 0; i < r2 -> L ; i++) {
    if(i%2) {
      ptr_dsread -> c2[j] |= dsbw1[(unsigned char)r2->line2[r2-> L - i -1]];
      j++; 
    } else {
      ptr_dsread -> c2[j] |= dsbw0[(unsigned char)r2->line2[r2 -> L - i -1]]; 
    }
  }
  // Look for the matches 
  uint64_t ul1, ul1sh, cmpl; 
  uint64_t ul2 = *((uint64_t *)ptr_dsread -> c2);
  int n; 
  int pos = r1 -> L - COMPACT*sizeof(uint64_t);
  int Nwindows = pos/2; 
  for (i = 0 ; i < Nwindows; i++) {
    ul1 = *(uint64_t *) (ptr_dsread -> c1 + Nwindows - i) ;
    ul1sh = *(uint64_t *) (ptr_dsread -> c1sh + Nwindows - i) ;
    cmpl = ~(ul1 ^ ul2); 
   n = __builtin_popcountll(cmpl);
    if (n > ptr_dsread -> nmax) {
      ptr_dsread -> nmax = n; 
      ptr_dsread -> posmax = pos;
    }
    pos--; 
    cmpl = ~(ul1sh ^ (ul2)); 
    n = __builtin_popcountll(cmpl);
    if (n > ptr_dsread -> nmax) {
      ptr_dsread -> nmax = n; 
      ptr_dsread -> posmax = pos;
    }
    pos--; 
  }
  
  // Check the 16 last nucleotides,
  if ( ptr_dsread -> nmax < 60) {
     uint32_t u1, u1sh, cmp;
     uint32_t u2 = *((uint32_t *)ptr_dsread -> c2);
     pos = r1 -> L - COMPACT*sizeof(uint32_t);
     Nwindows = pos/2; 
     for (i = 0 ; i < 8; i++) {
       u1 = *(uint32_t *) (ptr_dsread -> c1 + Nwindows - i) ;
       u1sh = *(uint32_t *) (ptr_dsread -> c1sh + Nwindows - i) ;
       cmp = ~(u1 ^ u2); 
       n = __builtin_popcountll(cmp);
       if (n > 28) {
         ptr_dsread -> nmax = 64;  //change afterwards
         ptr_dsread -> posmax = pos;
       }
       pos--; 
       cmp = ~(u1sh ^ (u2)); 
       n = __builtin_popcountll(cmp);
       if (n > 28) {
         ptr_dsread -> nmax = 64;  //change afterwards
         ptr_dsread -> posmax = pos;
       }
       pos--;
     }
  }
  // Compute score
}

