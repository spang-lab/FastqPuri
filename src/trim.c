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
 * @file trim.c
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 24.08.2017
 * @brief trims/filter sequences after Quality, N's contaminations.
 *
 * */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "trim.h"
#include "Lmer.h"
#include "defines.h"
#include "config.h"


extern char LT[256]; /**< global variable. Lookup table. */
extern char Nencode;
extern Iparam_trimFilter par_TF;

#define TRIM_STRING 20

/**
*
* @brief checks if a sequence contains any non standard base callings (N's)
* @returns 0 if no N's found, 1 if N's found
*
* This function checks if any of the base callings in a
* given fastq read is different from A, C, G, T. Basically,
* any char different from the former ones is classified as N.
* */
static int no_N(Fq_read *seq) {
  int i;
  char Lmer[seq -> L];
  memcpy(Lmer, seq -> line2, seq -> L);
  Lmer_sLmer(Lmer, seq -> L);
  for (i = 0; i < seq -> L; i++) {
    if (Lmer[i] >= Nencode) {
      return 0;
    }
  }
  return 1;
}

/**
 * @brief Finds the largest Nfree sub-seq and keeps it if larger than minL
 * @param seq fastq read
 * @param minL minimum accepted trimmed length
 * @return 0 if not used, 1 if accepted as is, 2 if accepted and trimmed
 *
 * */
static int Nfree_Lmer(Fq_read *seq, int minL) {
  int pos = 0;  // 1st pos. of the largest N-free sub-seq
  int pos_curr = 0;  // 1st pos. of the current N-free sub-seq (updated in loop)
  int len_max = 0;  // length of the largest N-free sub-seq
  int len_cur = 0;  // length of the current N-free sub-seq (updated in loop)
  int len_cum = 0;  // cumulative N-free length as we run in the loop
  int i;
  char Lmer[seq-> L];
  memcpy(Lmer, seq->line2, seq -> L);
  Lmer_sLmer(Lmer, seq->L);
  for (i = 0 ; i < seq -> L; i++) {
     len_cum++;
     if (Lmer[i] >= Nencode) {
        pos_curr = i + 1 - len_cum;
        len_cur = i - pos_curr;
        len_cum = 0;
        if (len_cur > len_max) {
           pos = pos_curr;
           len_max = len_cur;
        }
     }
  }
  if (len_cum == seq -> L) {
     return 1;
  }
  // Consider last piece
  if (len_cum > len_max) {
     len_max = len_cum;
     pos = seq -> L - len_max;
  }
  // Check the length, discard it if < Lmin, and trim it otherwise
  if (len_max < minL) {
     return 0;
  } else {
     memmove(seq -> line2, seq -> line2+pos, len_max);
     seq -> line2[len_max] = '\0';
     memmove(seq -> line4, seq -> line4+pos, len_max);
     seq -> line4[len_max] = '\0';
     seq -> L = len_max;
     char add[TRIM_STRING];
     snprintf(add, TRIM_STRING, " TRIMN:%d:%d", pos, pos+len_max);
     if (strlen(add) + strlen(seq -> line3) > READ_MAXLEN) {
       fprintf(stderr, "Cannot append  %s to %s.\n", add, seq -> line3);
       fprintf(stderr, "sequence exceeds the predifined limit.\n");
       fprintf(stderr, "Set READ_MAXLEN to a larger value in cmake\n");
       fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
       fprintf(stderr, "Exiting program.\n");
       exit(EXIT_FAILURE);
     }
     strncat(seq -> line3, add, TRIM_STRING);
     return 2;
  }
}

/**
 * @brief trims a read if N's are at the ends and the remaining sub-seq >= minL
 * @param seq fastq read
 * @param minL minimum accepted trimmed length
 * @return 0 if not used, 1 no N's found, 2 if accepted and trimmed
 *
 * */
static int Ntrim_ends(Fq_read *seq, int minL) {
  char Lmer[seq-> L];
  int t_start = 0;
  int t_end = seq -> L - 1;
  memcpy(Lmer, seq->line2, seq -> L);
  Lmer_sLmer(Lmer, seq->L);
  while (Lmer[t_start] >= Nencode)
     t_start++;
  while (Lmer[t_end] >= Nencode)
     t_end--;
  if ((t_end - t_start) == (seq -> L - 1)) {
     return 1;
  } else if ((t_end - t_start) < minL - 1) {
     return 0;
  } else {
    (seq -> L) = t_end - t_start + 1;
    memmove(seq -> line4, seq -> line4 + t_start, seq -> L);
    memmove(seq -> line2, seq -> line2 + t_start, seq -> L);
    seq -> line4[seq -> L] = '\0';
    seq -> line2[seq -> L] = '\0';
    char add[TRIM_STRING];
    snprintf(add, TRIM_STRING, " TRIMN:%d:%d", t_start, t_end);
    if (strlen(add) + strlen(seq -> line3) > READ_MAXLEN) {
      fprintf(stderr, "Cannot append  %s to %s.\n", add, seq -> line3);
      fprintf(stderr, "sequence exceeds the predifined limit.\n");
      fprintf(stderr, "Set READ_MAXLEN to a larger value in cmake\n");
      fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
      fprintf(stderr, "Exiting program.\n");
      exit(EXIT_FAILURE);
    }
    strncat(seq -> line3, add, TRIM_STRING);
    return 2;
  }
}

/**
 * @brief checks if a sequence contains lowQ nucleotides
 * @param seq fastq read
 * @param minQ minimum accepted quality value
 * @return 0 if seq contains lowQ nucleotides, 1 otherwise
 *
 * */
static int no_lowQ(Fq_read *seq, int minQ) {
  int i;
  for (i = 0; i < seq -> L; i++) {
    if (seq->line4[i] < (ZEROQ + minQ))
        return 0;
  }
  return 1;
}

/**
 * @brief trims a read if lowQs are at the ends and remaining sub-seq >= minL
 * @param minQ minimum accepted quality value
 * @param minL minimum accepted trimmed length
 * @return 0 if not used, 1 if accepted as is, 2 if accepted and trimmed
 *
 * */
static int Qtrim_ends(Fq_read *seq, int minQ, int minL) {
  // Beginning of sequence:
  int L = (seq->L)-1;
  int t_start = 0;
  int t_end = L;
  while (seq->line4[t_start] < (ZEROQ + minQ))
     t_start++;
  while (seq->line4[t_end] < (ZEROQ + minQ))
     t_end--;
  // Accept sequence as is
  if ((t_end - t_start) == L) {
     return 1;
  }
  // Discard sequence
  if ((t_end - t_start) < minL - 1) {
     return 0;
  }
  // Trim the sequence
  (seq -> L) = t_end - t_start + 1;
  memmove(seq -> line4, seq -> line4 + t_start, seq -> L);
  memmove(seq -> line2, seq -> line2 + t_start, seq -> L);
  seq -> line4[seq -> L] = '\0';
  seq -> line2[seq -> L] = '\0';
  char add[TRIM_STRING];
  snprintf(add, TRIM_STRING, " TRIMQ:%d:%d", t_start, t_end);
  if (strlen(add) + strlen(seq -> line3) > READ_MAXLEN) {
    fprintf(stderr, "Cannot append  %s to %s.\n", add, seq -> line3);
    fprintf(stderr, "sequence exceeds the predifined limit.\n");
    fprintf(stderr, "Set READ_MAXLEN to a larger value in cmake\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    fprintf(stderr, "Exiting program.\n");
    exit(EXIT_FAILURE);
  }
  strncat(seq -> line3, add, TRIM_STRING);
  return 2;
}

/**
 * @brief accepts the sequence as is if there are less than nlowQ
 * @param seq fastq read
 * @param minQ minimum accepted quality value
 * @param nlowQ threshold on lowQ nucleotides (>= NOT allowed)
 * @return 0 if not used, 1 if accepted as is
 * */
static int Qtrim_frac(Fq_read *seq, int minQ, int nlowQ) {
  int ilowQ = 0, i = 0;
  while (i < (seq->L) &&  ilowQ < nlowQ) {
     if (seq->line4[i] < (ZEROQ + minQ)) ilowQ++;
     i++;
  }
  return (ilowQ == nlowQ) ? 0: 1;
}

/**
 * @brief trims the ends for lowQ. The rest is kept if it contains < nlowQ lowQ
 * @param seq fastq read
 * @param minQ minimum accepted quality value
 * @param minL minimum accepted trimmed length
 * @param nlowQ threshold on lowQ nucleotides (>= NOT allowed)
 * @return 0 if not used, 1 if accepted as is
 *
 * Trims the ends if there are lowQ bases. From the remaining
 * part, it counts how many lowQ bases there are and keeps it if
 * there are less than nlowQ.
 *
 * */
static int Qtrim_endsfrac(Fq_read *seq, int minQ, int minL, int nlowQ ) {
  int L = (seq->L)-1;  // last accessible element of the sequence
  int t_start = 0;
  int t_end = L;
  int ilowQ = 0;
  int i;
  while (seq->line4[t_start] < (ZEROQ + minQ))
    t_start++;
  while (seq->line4[t_end] < (ZEROQ + minQ))
    t_end--;
  i = t_start;
  while (i < t_end && ilowQ < nlowQ) {
    if (seq->line4[i] < (ZEROQ + minQ)) ilowQ++;
    i++;
  }
  // Discard sequence
  if (ilowQ == nlowQ || (t_end - t_start) < minL - 1) return 0;
  // Accept sequence as is
  if ((t_end - t_start) == L) {
     return 1;
  }
  // Trim the sequence
  (seq -> L) = t_end - t_start + 1;
  memmove(seq -> line4, seq -> line4 + t_start, seq -> L);
  memmove(seq -> line2, seq -> line2 + t_start, seq -> L);
  seq -> line4[seq -> L] = '\0';
  seq -> line2[seq -> L] = '\0';
  char add[TRIM_STRING];
  snprintf(add, TRIM_STRING, " TRIMQ:%d:%d", t_start, t_end);
  if (strlen(add) + strlen(seq -> line3) > READ_MAXLEN) {
    fprintf(stderr, "Cannot append  %s to %s.\n", add, seq -> line3);
    fprintf(stderr, "sequence exceeds the predifined limit.\n");
    fprintf(stderr, "Set READ_MAXLEN to a larger value in cmake\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    fprintf(stderr, "Exiting program.\n");
    exit(EXIT_FAILURE);
  }
  strncat(seq -> line3, add, TRIM_STRING);
  return 2;
}

/**
 * @brief trims left from the left and right from the right
 * @param seq fastq read
 * @param left number of nucleotides to be trimmed from the left
 * @param right number of nucleotides to be trimmed from the right
 * @return 2, since they are all accepted and trim
 *
 * */
static int Qtrim_global(Fq_read *seq, int left, int right ) {
  int t_end = seq ->L - left;
  (seq -> L) -= (left+right);
  memmove(seq -> line4, seq -> line4 + left, seq -> L);
  memmove(seq -> line2, seq -> line2 + left, seq -> L);
  seq -> line4[seq -> L] = '\0';
  seq -> line2[seq -> L] = '\0';
  char add[TRIM_STRING];
  snprintf(add, TRIM_STRING, " TRIMQ:%d:%d", left, t_end);
  if (strlen(add) + strlen(seq -> line3) > READ_MAXLEN) {
    fprintf(stderr, "Cannot append  %s to %s.\n", add, seq -> line3);
    fprintf(stderr, "sequence exceeds the predifined limit.\n");
    fprintf(stderr, "Set READ_MAXLEN to a larger value in cmake\n");
    fprintf(stderr, "File: %s, line: %d\n", __FILE__, __LINE__);
    fprintf(stderr, "Exiting program.\n");
    exit(EXIT_FAILURE);
  }
  strncat(seq -> line3, add, TRIM_STRING);
  return 2;
}

/**
 * @brief trims sequence based on presence of N nucleotides
 * @param seq fastq read
 * @return -1 error, 0 discarded, 1 accepted as is, 2 accepted and trimmed
 *
 * This function calls a different function depending on the method
 * passed as input par_TF.trimN:
 * - NO(0):  accepts it as is, (1),
 * - ALL(1): accepts it as is if NO N's found (1), rejects it otherwise (0),
 * - ENDS(2): trims the ends and accepts it if it is longer than minL (2 if
 *            trimming, 1 if no trimming), rejects it otherwise (0),
 * - STRIP(3): finds the longest N-free subsequence and trims it if
 *             it is at least minL nucleotides long (2 if trimming, 1 if
 *             no N's are found), rejects it otherwise (0).
 * */
int trim_sequenceN(Fq_read *seq ) {
  return (par_TF.trimN == NO)? 1:
          (par_TF.trimN == ALL)? no_N(seq):
          (par_TF.trimN == ENDS)? Ntrim_ends(seq, par_TF.minL):
          (par_TF.trimN == STRIP)? Nfree_Lmer(seq, par_TF.minL): -1;
}

/**
 * @brief trims sequence based on lowQ base callings
 * @param seq fastq read
 * @return -1 error, 0 discarded, 1 accepted as is, 2 accepted and trimmed
 *
 * This function calls a different function depending on the method
 * passed as input par_TF.trimQ:
 * - NO(0): accepts is as is , (1),
 * - FRAC(1): accepts it if less than par_TF.nlowQ are found (1), rejects
 *            it otherwise (0),
 * - ENDS(2): trims the ends and accepts it if it is longer than minL  (2 if
 *            triming, 1 if no trimming), rejects it otherwise (0),
 * - ENDSFRAC(3): trims the ends and accepts if the remaining sequence
 *                is at least minL bases long and if it contains less than
 *                nlowQ lowQ nucleotides (2 if trimming, 1 if no trimming).
 *                Otherwise, it is rejected, (0).
 * - GLOBAL(4): it trims globally globleft nucleotides from the left and
 *              globright from the right, (returns 2).
 *
 * */
int trim_sequenceQ(Fq_read *seq) {
  return (par_TF.trimQ == NO)? 1 :
         (par_TF.trimQ == ALL)? no_lowQ(seq, par_TF.minQ):
         (par_TF.trimQ == ENDS) ? Qtrim_ends(seq, par_TF.minQ, par_TF.minL):
         (par_TF.trimQ == FRAC) ? Qtrim_frac(seq, par_TF.minQ, par_TF.nlowQ):
         (par_TF.trimQ == ENDSFRAC) ?
                Qtrim_endsfrac(seq, par_TF.minQ, par_TF.minL, par_TF.nlowQ):
         (par_TF.trimQ == GLOBAL) ?
                Qtrim_global(seq, par_TF.globleft, par_TF.globright): -1;
}

// Check if a read is in the sequence or its reverse complement
/**
 * @brief check if Lread is contained in tree.
 * @param tree_ptr pointer to Tree structure
 * @param seq fastq read
 * @param score threshold score [0,1] over which we consider a read was in tree
 * @returns true if read was found, false otherwise
 *
 * */
bool is_read_inTree(Tree *tree_ptr, Fq_read *seq) {
  char read[seq -> L];
  memcpy(read, seq -> line2, seq -> L);
  if (check_path(tree_ptr, read, seq -> L) > par_TF.score) {
     return true;
  } else {
     rev_comp(read, seq -> L);
     return (check_path(tree_ptr, read, seq -> L) > par_TF.score);
  }
}

