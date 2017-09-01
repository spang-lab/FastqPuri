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
 * @file str_manip.h
 * @author Paula Perez <paulaperezrubio@gmail.com>
 * @date 03.08.2017
 * @brief functions that do string manipulation
 *
 * */

#ifndef STR_MANIP_H_
#define STR_MANIP_H_


/**
 * @brief contains a splitted string and the number or splitted fields
 * */
typedef struct _split {
  int N;    /**< Number of substrings in which the string was splitted */ 
  char **s; /**< Substring array containing the splitted substrings */  
} Split; 

int strindex(char *s, char *t); 

int count_char(char *s, char c);

Split strsplit(char *str, char sep);

#endif  // endif STR_MANIP_H_
