/***************************************************************************
 *   Copyright (C) 2009 by Bharadwaj Subramanian   *
 *   bharadwajs@pupil.ices.utexas.edu   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or         *
 *   modify it under the terms of the GNU Lesser General Public            *
 *   License version 2.1 as published by the Free Software Foundation.     *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *   Lesser General Public License for more details.                       *
 *                                                                         *
 *   You should have received a copy of the GNU Lesser General Public      *
 *   License along with this library; if not, write to the Free Software   *
 *   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA         *
 *   02110-1301  USA                                                       *
 ***************************************************************************/
#include "MCTester.h"
#include "Grid.h"
#include "iostream"

using namespace std;


// gets the configuration that generated the case

    int confcase0[2][8] = {
  /* case  0 (  -1 ):   0 */   { 0, 0, 0, 0, 0, 0, 0, 0 },
  /* case  0 (  -1 ): 255 */   { 1, 1, 1, 1, 1, 1, 1, 1 }
    };

    int confcase1[16][8] = {
  /* case  1 (  0 ):   1 */   { 1, 0, 0, 0, 0, 0, 0, 0 },
  /* case  1 (  1 ):   2 */   { 0, 1, 0, 0, 0, 0, 0, 0 },
  /* case  1 (  2 ):   4 */   { 0, 0, 1, 0, 0, 0, 0, 0 },
  /* case  1 (  3 ):   8 */   { 0, 0, 0, 1, 0, 0, 0, 0 },
  /* case  1 (  4 ):  16 */   { 0, 0, 0, 0, 1, 0, 0, 0 },
  /* case  1 (  5 ):  32 */   { 0, 0, 0, 0, 0, 1, 0, 0 },
  /* case  1 (  6 ):  64 */   { 0, 0, 0, 0, 0, 0, 1, 0 },
  /* case  1 (  7 ): 128 */   { 0, 0, 0, 0, 0, 0, 0, 1 },
  /* case  1 (  8 ): 127 */   { 1, 1, 1, 1, 1, 1, 1, 0 },
  /* case  1 (  9 ): 191 */   { 1, 1, 1, 1, 1, 1, 0, 1 },
  /* case  1 ( 10 ): 223 */   { 1, 1, 1, 1, 1, 0, 1, 1 },
  /* case  1 ( 11 ): 239 */   { 1, 1, 1, 1, 0, 1, 1, 1 },
  /* case  1 ( 12 ): 247 */   { 1, 1, 1, 0, 1, 1, 1, 1 },
  /* case  1 ( 13 ): 251 */   { 1, 1, 0, 1, 1, 1, 1, 1 },
  /* case  1 ( 14 ): 253 */   { 1, 0, 1, 1, 1, 1, 1, 1 },
  /* case  1 ( 15 ): 254 */   { 0, 1, 1, 1, 1, 1, 1, 1 }
    };

    int confcase2[24][8] = {
  /* case  2 (  0 ):   3 */   { 1, 1, 0, 0, 0, 0, 0, 0 },
  /* case  2 (  1 ):   9 */   { 1, 0, 0, 1, 0, 0, 0, 0 },
  /* case  2 (  2 ):  17 */   { 1, 0, 0, 0, 1, 0, 0, 0 },
  /* case  2 (  3 ):   6 */   { 0, 1, 1, 0, 0, 0, 0, 0 },
  /* case  2 (  4 ):  34 */   { 0, 1, 0, 0, 0, 1, 0, 0 },
  /* case  2 (  5 ):  12 */   { 0, 0, 1, 1, 0, 0, 0, 0 },
  /* case  2 (  6 ):  68 */   { 0, 0, 1, 0, 0, 0, 1, 0 },
  /* case  2 (  7 ): 136 */   { 0, 0, 0, 1, 0, 0, 0, 1 },
  /* case  2 (  8 ):  48 */   { 0, 0, 0, 0, 1, 1, 0, 0 },
  /* case  2 (  9 ): 144 */   { 0, 0, 0, 0, 1, 0, 0, 1 },
  /* case  2 ( 10 ):  96 */   { 0, 0, 0, 0, 0, 1, 1, 0 },
  /* case  2 ( 11 ): 192 */   { 0, 0, 0, 0, 0, 0, 1, 1 },
  /* case  2 ( 12 ):  63 */   { 1, 1, 1, 1, 1, 1, 0, 0 },
  /* case  2 ( 13 ): 159 */   { 1, 1, 1, 1, 1, 0, 0, 1 },
  /* case  2 ( 14 ): 111 */   { 1, 1, 1, 1, 0, 1, 1, 0 },
  /* case  2 ( 15 ): 207 */   { 1, 1, 1, 1, 0, 0, 1, 1 },
  /* case  2 ( 16 ): 119 */   { 1, 1, 1, 0, 1, 1, 1, 0 },
  /* case  2 ( 17 ): 187 */   { 1, 1, 0, 1, 1, 1, 0, 1 },
  /* case  2 ( 18 ): 243 */   { 1, 1, 0, 0, 1, 1, 1, 1 },
  /* case  2 ( 19 ): 221 */   { 1, 0, 1, 1, 1, 0, 1, 1 },
  /* case  2 ( 20 ): 249 */   { 1, 0, 0, 1, 1, 1, 1, 1 },
  /* case  2 ( 21 ): 238 */   { 0, 1, 1, 1, 0, 1, 1, 1 },
  /* case  2 ( 22 ): 246 */   { 0, 1, 1, 0, 1, 1, 1, 1 },
  /* case  2 ( 23 ): 252 */   { 0, 0, 1, 1, 1, 1, 1, 1 }
    };

    int confcase3[24][8] = {
  /* case  3 (  0 ):   5 */   { 1, 0, 1, 0, 0, 0, 0, 0 },
  /* case  3 (  1 ):  33 */   { 1, 0, 0, 0, 0, 1, 0, 0 },
  /* case  3 (  2 ): 129 */   { 1, 0, 0, 0, 0, 0, 0, 1 },
  /* case  3 (  3 ):  10 */   { 0, 1, 0, 1, 0, 0, 0, 0 },
  /* case  3 (  4 ):  18 */   { 0, 1, 0, 0, 1, 0, 0, 0 },
  /* case  3 (  5 ):  66 */   { 0, 1, 0, 0, 0, 0, 1, 0 },
  /* case  3 (  6 ):  36 */   { 0, 0, 1, 0, 0, 1, 0, 0 },
  /* case  3 (  7 ): 132 */   { 0, 0, 1, 0, 0, 0, 0, 1 },
  /* case  3 (  8 ):  24 */   { 0, 0, 0, 1, 1, 0, 0, 0 },
  /* case  3 (  9 ):  72 */   { 0, 0, 0, 1, 0, 0, 1, 0 },
  /* case  3 ( 10 ):  80 */   { 0, 0, 0, 0, 1, 0, 1, 0 },
  /* case  3 ( 11 ): 160 */   { 0, 0, 0, 0, 0, 1, 0, 1 },
  /* case  3 ( 12 ):  95 */   { 1, 1, 1, 1, 1, 0, 1, 0 },
  /* case  3 ( 13 ): 175 */   { 1, 1, 1, 1, 0, 1, 0, 1 },
  /* case  3 ( 14 ): 183 */   { 1, 1, 1, 0, 1, 1, 0, 1 },
  /* case  3 ( 15 ): 231 */   { 1, 1, 1, 0, 0, 1, 1, 1 },
  /* case  3 ( 16 ): 123 */   { 1, 1, 0, 1, 1, 1, 1, 0 },
  /* case  3 ( 17 ): 219 */   { 1, 1, 0, 1, 1, 0, 1, 1 },
  /* case  3 ( 18 ): 189 */   { 1, 0, 1, 1, 1, 1, 0, 1 },
  /* case  3 ( 19 ): 237 */   { 1, 0, 1, 1, 0, 1, 1, 1 },
  /* case  3 ( 20 ): 245 */   { 1, 0, 1, 0, 1, 1, 1, 1 },
  /* case  3 ( 21 ): 126 */   { 0, 1, 1, 1, 1, 1, 1, 0 },
  /* case  3 ( 22 ): 222 */   { 0, 1, 1, 1, 1, 0, 1, 1 },
  /* case  3 ( 23 ): 250 */   { 0, 1, 0, 1, 1, 1, 1, 1 }
    };

    int confcase4[8][8] = {
  /* case  4 (  0 ):  65 */   { 1, 0, 0, 0, 0, 0, 1, 0 },
  /* case  4 (  1 ): 130 */   { 0, 1, 0, 0, 0, 0, 0, 1 },
  /* case  4 (  2 ):  20 */   { 0, 0, 1, 0, 1, 0, 0, 0 },
  /* case  4 (  3 ):  40 */   { 0, 0, 0, 1, 0, 1, 0, 0 },
  /* case  4 (  4 ): 215 */   { 1, 1, 1, 0, 1, 0, 1, 1 },
  /* case  4 (  5 ): 235 */   { 1, 1, 0, 1, 0, 1, 1, 1 },
  /* case  4 (  6 ): 125 */   { 1, 0, 1, 1, 1, 1, 1, 0 },
  /* case  4 (  7 ): 190 */   { 0, 1, 1, 1, 1, 1, 0, 1 }
    };

    int confcase5[48][8] = {
  /* case  5 (  0 ):   7 */   { 1, 1, 1, 0, 0, 0, 0, 0 },
  /* case  5 (  1 ):  11 */   { 1, 1, 0, 1, 0, 0, 0, 0 },
  /* case  5 (  2 ):  19 */   { 1, 1, 0, 0, 1, 0, 0, 0 },
  /* case  5 (  3 ):  35 */   { 1, 1, 0, 0, 0, 1, 0, 0 },
  /* case  5 (  4 ):  13 */   { 1, 0, 1, 1, 0, 0, 0, 0 },
  /* case  5 (  5 ):  25 */   { 1, 0, 0, 1, 1, 0, 0, 0 },
  /* case  5 (  6 ): 137 */   { 1, 0, 0, 1, 0, 0, 0, 1 },
  /* case  5 (  7 ):  49 */   { 1, 0, 0, 0, 1, 1, 0, 0 },
  /* case  5 (  8 ): 145 */   { 1, 0, 0, 0, 1, 0, 0, 1 },
  /* case  5 (  9 ):  14 */   { 0, 1, 1, 1, 0, 0, 0, 0 },
  /* case  5 ( 10 ):  38 */   { 0, 1, 1, 0, 0, 1, 0, 0 },
  /* case  5 ( 11 ):  70 */   { 0, 1, 1, 0, 0, 0, 1, 0 },
  /* case  5 ( 12 ):  50 */   { 0, 1, 0, 0, 1, 1, 0, 0 },
  /* case  5 ( 13 ):  98 */   { 0, 1, 0, 0, 0, 1, 1, 0 },
  /* case  5 ( 14 ):  76 */   { 0, 0, 1, 1, 0, 0, 1, 0 },
  /* case  5 ( 15 ): 140 */   { 0, 0, 1, 1, 0, 0, 0, 1 },
  /* case  5 ( 16 ): 100 */   { 0, 0, 1, 0, 0, 1, 1, 0 },
  /* case  5 ( 17 ): 196 */   { 0, 0, 1, 0, 0, 0, 1, 1 },
  /* case  5 ( 18 ): 152 */   { 0, 0, 0, 1, 1, 0, 0, 1 },
  /* case  5 ( 19 ): 200 */   { 0, 0, 0, 1, 0, 0, 1, 1 },
  /* case  5 ( 20 ): 112 */   { 0, 0, 0, 0, 1, 1, 1, 0 },
  /* case  5 ( 21 ): 176 */   { 0, 0, 0, 0, 1, 1, 0, 1 },
  /* case  5 ( 22 ): 208 */   { 0, 0, 0, 0, 1, 0, 1, 1 },
  /* case  5 ( 23 ): 224 */   { 0, 0, 0, 0, 0, 1, 1, 1 },
  /* case  5 ( 24 ):  31 */   { 1, 1, 1, 1, 1, 0, 0, 0 },
  /* case  5 ( 25 ):  47 */   { 1, 1, 1, 1, 0, 1, 0, 0 },
  /* case  5 ( 26 ):  79 */   { 1, 1, 1, 1, 0, 0, 1, 0 },
  /* case  5 ( 27 ): 143 */   { 1, 1, 1, 1, 0, 0, 0, 1 },
  /* case  5 ( 28 ):  55 */   { 1, 1, 1, 0, 1, 1, 0, 0 },
  /* case  5 ( 29 ): 103 */   { 1, 1, 1, 0, 0, 1, 1, 0 },
  /* case  5 ( 30 ):  59 */   { 1, 1, 0, 1, 1, 1, 0, 0 },
  /* case  5 ( 31 ): 155 */   { 1, 1, 0, 1, 1, 0, 0, 1 },
  /* case  5 ( 32 ): 115 */   { 1, 1, 0, 0, 1, 1, 1, 0 },
  /* case  5 ( 33 ): 179 */   { 1, 1, 0, 0, 1, 1, 0, 1 },
  /* case  5 ( 34 ): 157 */   { 1, 0, 1, 1, 1, 0, 0, 1 },
  /* case  5 ( 35 ): 205 */   { 1, 0, 1, 1, 0, 0, 1, 1 },
  /* case  5 ( 36 ): 185 */   { 1, 0, 0, 1, 1, 1, 0, 1 },
  /* case  5 ( 37 ): 217 */   { 1, 0, 0, 1, 1, 0, 1, 1 },
  /* case  5 ( 38 ): 241 */   { 1, 0, 0, 0, 1, 1, 1, 1 },
  /* case  5 ( 39 ): 110 */   { 0, 1, 1, 1, 0, 1, 1, 0 },
  /* case  5 ( 40 ): 206 */   { 0, 1, 1, 1, 0, 0, 1, 1 },
  /* case  5 ( 41 ): 118 */   { 0, 1, 1, 0, 1, 1, 1, 0 },
  /* case  5 ( 42 ): 230 */   { 0, 1, 1, 0, 0, 1, 1, 1 },
  /* case  5 ( 43 ): 242 */   { 0, 1, 0, 0, 1, 1, 1, 1 },
  /* case  5 ( 44 ): 220 */   { 0, 0, 1, 1, 1, 0, 1, 1 },
  /* case  5 ( 45 ): 236 */   { 0, 0, 1, 1, 0, 1, 1, 1 },
  /* case  5 ( 46 ): 244 */   { 0, 0, 1, 0, 1, 1, 1, 1 },
  /* case  5 ( 47 ): 248 */   { 0, 0, 0, 1, 1, 1, 1, 1 }
    };

    int confcase6[48][8] = {
  /* case  6 (  0 ):  67 */   { 1, 1, 0, 0, 0, 0, 1, 0 },
  /* case  6 (  1 ): 131 */   { 1, 1, 0, 0, 0, 0, 0, 1 },
  /* case  6 (  2 ):  21 */   { 1, 0, 1, 0, 1, 0, 0, 0 },
  /* case  6 (  3 ):  69 */   { 1, 0, 1, 0, 0, 0, 1, 0 },
  /* case  6 (  4 ):  41 */   { 1, 0, 0, 1, 0, 1, 0, 0 },
  /* case  6 (  5 ):  73 */   { 1, 0, 0, 1, 0, 0, 1, 0 },
  /* case  6 (  6 ):  81 */   { 1, 0, 0, 0, 1, 0, 1, 0 },
  /* case  6 (  7 ):  97 */   { 1, 0, 0, 0, 0, 1, 1, 0 },
  /* case  6 (  8 ): 193 */   { 1, 0, 0, 0, 0, 0, 1, 1 },
  /* case  6 (  9 ):  22 */   { 0, 1, 1, 0, 1, 0, 0, 0 },
  /* case  6 ( 10 ): 134 */   { 0, 1, 1, 0, 0, 0, 0, 1 },
  /* case  6 ( 11 ):  42 */   { 0, 1, 0, 1, 0, 1, 0, 0 },
  /* case  6 ( 12 ): 138 */   { 0, 1, 0, 1, 0, 0, 0, 1 },
  /* case  6 ( 13 ): 146 */   { 0, 1, 0, 0, 1, 0, 0, 1 },
  /* case  6 ( 14 ): 162 */   { 0, 1, 0, 0, 0, 1, 0, 1 },
  /* case  6 ( 15 ): 194 */   { 0, 1, 0, 0, 0, 0, 1, 1 },
  /* case  6 ( 16 ):  28 */   { 0, 0, 1, 1, 1, 0, 0, 0 },
  /* case  6 ( 17 ):  44 */   { 0, 0, 1, 1, 0, 1, 0, 0 },
  /* case  6 ( 18 ):  52 */   { 0, 0, 1, 0, 1, 1, 0, 0 },
  /* case  6 ( 19 ):  84 */   { 0, 0, 1, 0, 1, 0, 1, 0 },
  /* case  6 ( 20 ): 148 */   { 0, 0, 1, 0, 1, 0, 0, 1 },
  /* case  6 ( 21 ):  56 */   { 0, 0, 0, 1, 1, 1, 0, 0 },
  /* case  6 ( 22 ): 104 */   { 0, 0, 0, 1, 0, 1, 1, 0 },
  /* case  6 ( 23 ): 168 */   { 0, 0, 0, 1, 0, 1, 0, 1 },
  /* case  6 ( 24 ):  87 */   { 1, 1, 1, 0, 1, 0, 1, 0 },
  /* case  6 ( 25 ): 151 */   { 1, 1, 1, 0, 1, 0, 0, 1 },
  /* case  6 ( 26 ): 199 */   { 1, 1, 1, 0, 0, 0, 1, 1 },
  /* case  6 ( 27 ): 107 */   { 1, 1, 0, 1, 0, 1, 1, 0 },
  /* case  6 ( 28 ): 171 */   { 1, 1, 0, 1, 0, 1, 0, 1 },
  /* case  6 ( 29 ): 203 */   { 1, 1, 0, 1, 0, 0, 1, 1 },
  /* case  6 ( 30 ): 211 */   { 1, 1, 0, 0, 1, 0, 1, 1 },
  /* case  6 ( 31 ): 227 */   { 1, 1, 0, 0, 0, 1, 1, 1 },
  /* case  6 ( 32 ):  61 */   { 1, 0, 1, 1, 1, 1, 0, 0 },
  /* case  6 ( 33 ):  93 */   { 1, 0, 1, 1, 1, 0, 1, 0 },
  /* case  6 ( 34 ): 109 */   { 1, 0, 1, 1, 0, 1, 1, 0 },
  /* case  6 ( 35 ): 117 */   { 1, 0, 1, 0, 1, 1, 1, 0 },
  /* case  6 ( 36 ): 213 */   { 1, 0, 1, 0, 1, 0, 1, 1 },
  /* case  6 ( 37 ): 121 */   { 1, 0, 0, 1, 1, 1, 1, 0 },
  /* case  6 ( 38 ): 233 */   { 1, 0, 0, 1, 0, 1, 1, 1 },
  /* case  6 ( 39 ):  62 */   { 0, 1, 1, 1, 1, 1, 0, 0 },
  /* case  6 ( 40 ): 158 */   { 0, 1, 1, 1, 1, 0, 0, 1 },
  /* case  6 ( 41 ): 174 */   { 0, 1, 1, 1, 0, 1, 0, 1 },
  /* case  6 ( 42 ): 182 */   { 0, 1, 1, 0, 1, 1, 0, 1 },
  /* case  6 ( 43 ): 214 */   { 0, 1, 1, 0, 1, 0, 1, 1 },
  /* case  6 ( 44 ): 186 */   { 0, 1, 0, 1, 1, 1, 0, 1 },
  /* case  6 ( 45 ): 234 */   { 0, 1, 0, 1, 0, 1, 1, 1 },
  /* case  6 ( 46 ): 124 */   { 0, 0, 1, 1, 1, 1, 1, 0 },
  /* case  6 ( 47 ): 188 */   { 0, 0, 1, 1, 1, 1, 0, 1 }
    };

    int confcase7[16][8] = {
  /* case  7 (  0 ):  37 */   { 1, 0, 1, 0, 0, 1, 0, 0 },
  /* case  7 (  1 ): 133 */   { 1, 0, 1, 0, 0, 0, 0, 1 },
  /* case  7 (  2 ): 161 */   { 1, 0, 0, 0, 0, 1, 0, 1 },
  /* case  7 (  3 ):  26 */   { 0, 1, 0, 1, 1, 0, 0, 0 },
  /* case  7 (  4 ):  74 */   { 0, 1, 0, 1, 0, 0, 1, 0 },
  /* case  7 (  5 ):  82 */   { 0, 1, 0, 0, 1, 0, 1, 0 },
  /* case  7 (  6 ): 164 */   { 0, 0, 1, 0, 0, 1, 0, 1 },
  /* case  7 (  7 ):  88 */   { 0, 0, 0, 1, 1, 0, 1, 0 },
  /* case  7 (  8 ): 167 */   { 1, 1, 1, 0, 0, 1, 0, 1 },
  /* case  7 (  9 ):  91 */   { 1, 1, 0, 1, 1, 0, 1, 0 },
  /* case  7 ( 10 ): 173 */   { 1, 0, 1, 1, 0, 1, 0, 1 },
  /* case  7 ( 11 ): 181 */   { 1, 0, 1, 0, 1, 1, 0, 1 },
  /* case  7 ( 12 ): 229 */   { 1, 0, 1, 0, 0, 1, 1, 1 },
  /* case  7 ( 13 ):  94 */   { 0, 1, 1, 1, 1, 0, 1, 0 },
  /* case  7 ( 14 ): 122 */   { 0, 1, 0, 1, 1, 1, 1, 0 },
  /* case  7 ( 15 ): 218 */   { 0, 1, 0, 1, 1, 0, 1, 1 }
    };

    int confcase8[6][8] = {
  /* case  8 (  0 ):  15 */   { 1, 1, 1, 1, 0, 0, 0, 0 },
  /* case  8 (  1 ):  51 */   { 1, 1, 0, 0, 1, 1, 0, 0 },
  /* case  8 (  2 ): 153 */   { 1, 0, 0, 1, 1, 0, 0, 1 },
  /* case  8 (  3 ): 102 */   { 0, 1, 1, 0, 0, 1, 1, 0 },
  /* case  8 (  4 ): 204 */   { 0, 0, 1, 1, 0, 0, 1, 1 },
  /* case  8 (  5 ): 240 */   { 0, 0, 0, 0, 1, 1, 1, 1 }
    };

    int confcase9[8][8] = {
  /* case  9 (  0 ):  39 */   { 1, 1, 1, 0, 0, 1, 0, 0 },
  /* case  9 (  1 ):  27 */   { 1, 1, 0, 1, 1, 0, 0, 0 },
  /* case  9 (  2 ): 141 */   { 1, 0, 1, 1, 0, 0, 0, 1 },
  /* case  9 (  3 ): 177 */   { 1, 0, 0, 0, 1, 1, 0, 1 },
  /* case  9 (  4 ):  78 */   { 0, 1, 1, 1, 0, 0, 1, 0 },
  /* case  9 (  5 ): 114 */   { 0, 1, 0, 0, 1, 1, 1, 0 },
  /* case  9 (  6 ): 228 */   { 0, 0, 1, 0, 0, 1, 1, 1 },
  /* case  9 (  7 ): 216 */   { 0, 0, 0, 1, 1, 0, 1, 1 }
    };

    int confcase10[6][8] = {
  /* case 10 (  0 ): 195 */   { 1, 1, 0, 0, 0, 0, 1, 1 },
  /* case 10 (  1 ):  85 */   { 1, 0, 1, 0, 1, 0, 1, 0 },
  /* case 10 (  2 ): 105 */   { 1, 0, 0, 1, 0, 1, 1, 0 },
  /* case 10 (  3 ): 150 */   { 0, 1, 1, 0, 1, 0, 0, 1 },
  /* case 10 (  4 ): 170 */   { 0, 1, 0, 1, 0, 1, 0, 1 },
  /* case 10 (  5 ):  60 */   { 0, 0, 1, 1, 1, 1, 0, 0 }
    };

    int confcase11[12][8] = {
  /* case 11 (  0 ):  23 */   { 1, 1, 1, 0, 1, 0, 0, 0 },
  /* case 11 (  1 ): 139 */   { 1, 1, 0, 1, 0, 0, 0, 1 },
  /* case 11 (  2 ):  99 */   { 1, 1, 0, 0, 0, 1, 1, 0 },
  /* case 11 (  3 ):  77 */   { 1, 0, 1, 1, 0, 0, 1, 0 },
  /* case 11 (  4 ):  57 */   { 1, 0, 0, 1, 1, 1, 0, 0 },
  /* case 11 (  5 ): 209 */   { 1, 0, 0, 0, 1, 0, 1, 1 },
  /* case 11 (  6 ):  46 */   { 0, 1, 1, 1, 0, 1, 0, 0 },
  /* case 11 (  7 ): 198 */   { 0, 1, 1, 0, 0, 0, 1, 1 },
  /* case 11 (  8 ): 178 */   { 0, 1, 0, 0, 1, 1, 0, 1 },
  /* case 11 (  9 ): 156 */   { 0, 0, 1, 1, 1, 0, 0, 1 },
  /* case 11 ( 10 ): 116 */   { 0, 0, 1, 0, 1, 1, 1, 0 },
  /* case 11 ( 11 ): 232 */   { 0, 0, 0, 1, 0, 1, 1, 1 }
    };

    int confcase12[24][8] = {
  /* case 12 (  0 ): 135 */   { 1, 1, 1, 0, 0, 0, 0, 1 },
  /* case 12 (  1 ):  75 */   { 1, 1, 0, 1, 0, 0, 1, 0 },
  /* case 12 (  2 ):  83 */   { 1, 1, 0, 0, 1, 0, 1, 0 },
  /* case 12 (  3 ): 163 */   { 1, 1, 0, 0, 0, 1, 0, 1 },
  /* case 12 (  4 ):  45 */   { 1, 0, 1, 1, 0, 1, 0, 0 },
  /* case 12 (  5 ):  53 */   { 1, 0, 1, 0, 1, 1, 0, 0 },
  /* case 12 (  6 ): 149 */   { 1, 0, 1, 0, 1, 0, 0, 1 },
  /* case 12 (  7 ): 101 */   { 1, 0, 1, 0, 0, 1, 1, 0 },
  /* case 12 (  8 ): 197 */   { 1, 0, 1, 0, 0, 0, 1, 1 },
  /* case 12 (  9 ):  89 */   { 1, 0, 0, 1, 1, 0, 1, 0 },
  /* case 12 ( 10 ): 169 */   { 1, 0, 0, 1, 0, 1, 0, 1 },
  /* case 12 ( 11 ): 225 */   { 1, 0, 0, 0, 0, 1, 1, 1 },
  /* case 12 ( 12 ):  30 */   { 0, 1, 1, 1, 1, 0, 0, 0 },
  /* case 12 ( 13 ):  86 */   { 0, 1, 1, 0, 1, 0, 1, 0 },
  /* case 12 ( 14 ): 166 */   { 0, 1, 1, 0, 0, 1, 0, 1 },
  /* case 12 ( 15 ):  58 */   { 0, 1, 0, 1, 1, 1, 0, 0 },
  /* case 12 ( 16 ): 154 */   { 0, 1, 0, 1, 1, 0, 0, 1 },
  /* case 12 ( 17 ): 106 */   { 0, 1, 0, 1, 0, 1, 1, 0 },
  /* case 12 ( 18 ): 202 */   { 0, 1, 0, 1, 0, 0, 1, 1 },
  /* case 12 ( 19 ): 210 */   { 0, 1, 0, 0, 1, 0, 1, 1 },
  /* case 12 ( 20 ):  92 */   { 0, 0, 1, 1, 1, 0, 1, 0 },
  /* case 12 ( 21 ): 172 */   { 0, 0, 1, 1, 0, 1, 0, 1 },
  /* case 12 ( 22 ): 180 */   { 0, 0, 1, 0, 1, 1, 0, 1 },
  /* case 12 ( 23 ): 120 */   { 0, 0, 0, 1, 1, 1, 1, 0 }
    };

    int confcase13[2][8] = {
  /* case 13 (  0 ): 165 */   { 1, 0, 1, 0, 0, 1, 0, 1 },
  /* case 13 (  1 ):  90 */   { 0, 1, 0, 1, 1, 0, 1, 0 }
    };

    int confcase14[12][8] = {
  /* case 14 (  0 ):  71 */   { 1, 1, 1, 0, 0, 0, 1, 0 },
  /* case 14 (  1 ):  43 */   { 1, 1, 0, 1, 0, 1, 0, 0 },
  /* case 14 (  2 ): 147 */   { 1, 1, 0, 0, 1, 0, 0, 1 },
  /* case 14 (  3 ):  29 */   { 1, 0, 1, 1, 1, 0, 0, 0 },
  /* case 14 (  4 ): 201 */   { 1, 0, 0, 1, 0, 0, 1, 1 },
  /* case 14 (  5 ): 113 */   { 1, 0, 0, 0, 1, 1, 1, 0 },
  /* case 14 (  6 ): 142 */   { 0, 1, 1, 1, 0, 0, 0, 1 },
  /* case 14 (  7 ):  54 */   { 0, 1, 1, 0, 1, 1, 0, 0 },
  /* case 14 (  8 ): 226 */   { 0, 1, 0, 0, 0, 1, 1, 1 },
  /* case 14 (  9 ): 108 */   { 0, 0, 1, 1, 0, 1, 1, 0 },
  /* case 14 ( 10 ): 212 */   { 0, 0, 1, 0, 1, 0, 1, 1 },
  /* case 14 ( 11 ): 184 */   { 0, 0, 0, 1, 1, 1, 0, 1 }
    };
//________________________________________________

    const int confsizes[15] = {
      2, 16, 24, 24, 8, 48, 48, 16, 6, 8, 6, 12, 24, 2, 12
    };
    
int checkFaceSaddle(float *val,int face_idx,float iso_val)
{
  float f[4]; // Variable holding the face vertex values.
  char v_signs[4]; // The signs of the vertices on the face.
  float fsaddle;

  f[0]=val[(int)faces[face_idx][0]];
  f[1]=val[(int)faces[face_idx][1]];
  f[2]=val[(int)faces[face_idx][2]];
  f[3]=val[(int)faces[face_idx][3]];

  v_signs[0]=(f[0]>iso_val?0:1);
  v_signs[1]=(f[1]>iso_val?0:1);
  v_signs[2]=(f[2]>iso_val?0:1);
  v_signs[3]=(f[3]>iso_val?0:1);
  
  // Now that we've got the right face values populated, we first check if
  // two opposing vertices have same signs.
  if(v_signs[0]==v_signs[2]&&v_signs[1]==v_signs[3]&&v_signs[0]!=v_signs[1])
  {
    fsaddle=(f[0]*f[2]-f[1]*f[3])/(f[0]-f[1]+f[2]-f[3]);
    // Now check if fsaddle>iso or <iso.

    if(fsaddle>iso_val) return 0; // Outside
    else /*if (fsaddle<iso_val)*/ return 1; // Inside
//     else
//     {
//       // 0 case; return the sign of a-b+c-d.
//       return (f[0]-f[1]+f[2]-f[3]<0);
//     }
  }
  else
  {
    return -1; // Don't know, don't care.
  }
}

bool isInCube(float* temp_vtx)
{
  int i;
  for (i=0;i<3;i++) {
    if (temp_vtx[i]<=0 || temp_vtx[i]>=1) return false;
  }
  return true;
}

float getTriVal(float val[8], float x, float y, float z, int res)
{
  float x_ratio,y_ratio,z_ratio;
  float temp1,temp2,temp3,temp4,temp5,temp6;

  x_ratio=((float)(res-x))/((float)res);
  y_ratio=((float)(res-y))/((float)res);
  z_ratio=((float)(res-z))/((float)res);
            
    
  temp1 = val[1] + (val[0]-val[1])*x_ratio;
  temp2 = val[3] + (val[2]-val[3])*x_ratio;
  temp3 = val[5] + (val[4]-val[5])*x_ratio;
  temp4 = val[7] + (val[6]-val[7])*x_ratio;
    

  temp5 = temp2  + (temp1-temp2)*y_ratio;
  temp6 = temp4  + (temp3-temp4)*y_ratio;

  return temp6  + (temp5-temp6)*z_ratio;
}

int computeBodySaddle(float* bval,float iso_val,vector<char> &bsaddlesigns)
{
  int i;
  float a,b,c,d,e,f,g,h,k0,k1,k2,z1,z2;
  float A,B,C,D,E,F,G,H;
  float bsaddle1[3],bsaddle2[3];


  A=bval[0];
  B=bval[1];
  C=bval[2];
  D=bval[3];
  E=bval[4];
  F=bval[5];
  G=bval[6];
  H=bval[7];

  a=A;
  b=-A+E;
  c=-A+D;
  d=A-D-E+H;
  e=-A+B;
  f=A-B-E+F;
  g=A-B+C-D;
  h=-A+B-C+D+E-F+G-H;

  /*
  h= f000;
  e=-f000+f100;
  f=-f000+f001;
  g=-f000+f110;
  b=f000-f001-f100+f101;
  d=f000-f110-f100+f010;
  c=f000-f110-f001+f111;
  a=-f000+f110+f001-f111
      +f100-f010-f101+f011;
 */


//   k0=c*f-b*g;
//   k1=d*f-b*h;
//   k2=d*g-c*h;
// 
//   float temp1,temp2,temp3;
//   temp1=g*g*k1*k1;
//   temp2=h*sqrt(k1)*(e*k2+g*k0);
//   temp3=temp1-temp2;
// 
//   if(temp3<0)
//   {
//     // We have a problem!
//     float bsaddleval=0;
//     for(int i=0;i<8;i++) bsaddleval+=bval[i];
//     bsaddleval/=8;
//     if(bsaddleval>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
// 
//     return 1;
//   }
// 
//   
//   float zsqrt=sqrt(temp3);
// 
//   z1=(zsqrt-g)/h;
//   z2=-1*(zsqrt+g)/h;
// 
//   bsaddle1[0]= -1*(c+d*z1)/(g+h*z1);
//   bsaddle1[1]= (k0+k1*z1)/k2;
//   bsaddle1[2]= z1;
// 
//   bsaddle2[0]= -1*(c+d*z2)/(g+h*z2);
//   bsaddle2[1]= (k0+k1*z2)/k2;
//   bsaddle2[2]= z2;

  k0=h*f*d-h*h*b;
  k1=2*g*(f*d-b*h);
  k2=g*f*c-g*g*b-e*h*c+e*d*g;

  if((k1*k1-4*k0*k2)<0)
  {
    float bsaddleval=0;
    for(int i=0;i<8;i++) bsaddleval+=bval[i];
    bsaddleval/=8;
    if(bsaddleval>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);

    cout<<"Imaginary body saddle found. Using center as body saddle. "<<endl;

    return 1;
  }

  float zsqrt=sqrt(k1*k1-4*k0*k2);
  z1=(-1*k1+zsqrt)/(2*k0);
  z2=(-1*k1-zsqrt)/(2*k0);

  bsaddle1[0]=-(d*z1+c)/(h*z1+g);
  bsaddle1[1]=(f*d*z1+f*c-b*g-b*h*z1)/(d*g-h*c);
  bsaddle1[2]=z1;

  bsaddle2[0]=-(d*z2+c)/(h*z2+g);
  bsaddle2[1]=(f*d*z2+f*c-b*g-b*h*z2)/(d*g-h*c);
  bsaddle2[2]=z2;
  
  float temp_vtx[3];


  if (isInCube(bsaddle1)) {
    if (isInCube(bsaddle2)) {
      //getTriVal() interpolates using trilinear interpolation.
      if (getTriVal(bval,bsaddle1[0],bsaddle1[1],bsaddle1[2],1) < getTriVal(bval,bsaddle2[0],bsaddle2[1],bsaddle2[2],1)) {
        for (i=0;i<3;i++) temp_vtx[i]=bsaddle2[i];
        for (i=0;i<3;i++) bsaddle2[i]=bsaddle1[i];
        for (i=0;i<3;i++) bsaddle1[i]=temp_vtx[i];
      }
      if(getTriVal(bval,bsaddle1[0],bsaddle1[1],bsaddle1[2],1)>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
      if(getTriVal(bval,bsaddle2[0],bsaddle2[1],bsaddle2[2],1)>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
      return 2;
    }
    if(getTriVal(bval,bsaddle1[0],bsaddle1[1],bsaddle1[2],1)>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
    return 1;
  }
  else {
    if (isInCube(bsaddle2)) {
      for (i=0;i<3;i++) temp_vtx[i]=bsaddle2[i];
      for (i=0;i<3;i++) bsaddle2[i]=bsaddle1[i];
      for (i=0;i<3;i++) bsaddle1[i]=temp_vtx[i];
      if(getTriVal(bval,bsaddle1[0],bsaddle1[1],bsaddle1[2],1)>=iso_val) bsaddlesigns.push_back(0); else bsaddlesigns.push_back(1);
      return 1;
    }
    return 0;
  }

}
//*
bool checkCases(int *given, float *vals, float iso_val, int casesSize)
{
  for(int i=0;i<casesSize;i++)
  {
    int j=0;
    // Test vertex signs and face saddles.
    for(;j<14;j++)
      if(mcCases[i][j]!=given[j]) break;
    if(j==14)
    {
      // See if body saddle needs to be checked.
      if(mcCases[i][14]!=-1)
      {
        // Check if required.
        
        vector<char> bsaddlesigns;
        int temp[8]={0,1,0,1,0,0,0,1};
        bool res=true;
        for(int t=0;t<8;t++) if(temp[t]!=given[t]) res=false;

        if(res)
        {
          int b=0;
        }
        computeBodySaddle(vals,iso_val,bsaddlesigns);
        // If no bodysaddles are found (we need one!) we fix by using the center.
        if(bsaddlesigns.size()==0)
        {
          for(int t=0;t<8;t++)
            cout<<vals[t]<<",";

          for(int t=0;t<14;t++)
            cout<<given[t]<<",";

          cout<<"No body saddles found"<<endl;

            float bsaddleval=0;
            for(int i=0;i<8;i++) bsaddleval+=vals[i];
            bsaddleval/=8;
            if(bsaddleval>=iso_val) given[14]=0; else given[14]=1;
        }
        else if(bsaddlesigns.size()==1)
          given[14]=bsaddlesigns[0];
        else if(bsaddlesigns.size()==2)
          if(bsaddlesigns[0]==bsaddlesigns[1])
            given[14]=0;
          else
           given[14]=1;
        
        if(mcCases[i][14]==given[14])
        {
            for(j=15;j<21;j++)
                given[j]=mcCases[i][j];
            return true;
        }
      }
      else
      {
        for(j=15;j<21;j++)
          given[j]=mcCases[i][j];
        return true;
      }
    }
  }
  return false;
}

MCCase MCTester::identifyMCCase(Cell &values,float iso_val)
{
 float *vals=values.vertexValues;
  int signs[21];
  for(int i=0;i<15;i++) signs[i]=-1;

  // Set the vertex signs
  for(int i=0;i<8;i++) if(vals[i]>=iso_val) signs[i] = 0; else signs[i] = 1;

  // Set the face saddle signs
  for(int i=8;i<14;i++)
  {
      int j=i-8;
      if( signs[faces[j][0]]==signs[faces[j][2]]&&
          signs[faces[j][1]]==signs[faces[j][3]]&&
          signs[faces[j][0]]!=signs[faces[j][1]])
      {
        signs[i]=checkFaceSaddle(vals,j,iso_val);
      }
  }
  bool res=true;
  int arr[8]={0,1,1,1,1,0,0,1};

  for(int t=0;t<8;t++)
    if(arr[t]!=signs[t]) res=false;

  if(res)
  {
    int b=0;
  }
  
  // Set the body saddle signs
  vector<char> bsaddlesigns;
//   computeBodySaddle(vals,iso_val,bsaddlesigns);
//   if(bsaddlesigns.size()==1)
//     signs[14]=bsaddlesigns[0];
//   else if(bsaddlesigns.size()==2)
//     if(bsaddlesigns[0]==bsaddlesigns[1])
//       signs[14]=0;
//     else
//       signs[14]=1;

    if(bsaddlesigns.size()==0)
    {
      int a=0;
    }
//   else
//   {
    // We *require* a body saddle to be present. So we check the center's sign.
//     float bsaddleval=0;
//     for(int i=0;i<8;i++) bsaddleval+=vals[i];
//     bsaddleval/=8;
//     if(bsaddleval>=iso_val) signs[14]=0; else signs[14]=1;
//   }

  bool result=true;

  
  result=checkCases(signs,vals,iso_val, 730);

  // Populate the right data structure.
  MCCase mcCase;
  mcCase.mcCase=signs[15];
  mcCase.faceIndex=signs[16];
  mcCase.bodyIndex=signs[17];
  mcCase.caseIndex=signs[20];

  mcCase.componentEdges=getComponentEdges(signs);

  return mcCase; 

}
/*/

bool checkCases(int cases[50][8],int given[8],int casesSize)
{
  for(int i=0;i<casesSize;i++)
  {
    int j=0;
    for(;j<8;j++)
      if(cases[i][j]!=given[j]) break;
    if(j==8)
      return true;
  }
  return false;
}

MCCase MCTester::identifyMCCase(Cell &values,float iso_val)
{
  float *vals=values.vertexValues;
  int oc_signs[8];
  int confcase=0;
  bool caseFound=false;

  MCCase mccase;

  for(int i=0;i<8;i++) if(vals[i]>iso_val) oc_signs[i] = 0; else oc_signs[i] = 1;
  

  if(!caseFound&&checkCases(confcase0,oc_signs,confsizes[0])) { caseFound=true; confcase=0; }
  if(!caseFound&&checkCases(confcase1,oc_signs,confsizes[1])) { caseFound=true; confcase=1; }
  if(!caseFound&&checkCases(confcase2,oc_signs,confsizes[2])) { caseFound=true; confcase=2; }
  if(!caseFound&&checkCases(confcase3,oc_signs,confsizes[3])) { caseFound=true; confcase=3; }
  if(!caseFound&&checkCases(confcase4,oc_signs,confsizes[4])) { caseFound=true; confcase=4; }
  if(!caseFound&&checkCases(confcase5,oc_signs,confsizes[5])) { caseFound=true; confcase=5; }
  if(!caseFound&&checkCases(confcase6,oc_signs,confsizes[6])) { caseFound=true; confcase=6; }
  if(!caseFound&&checkCases(confcase7,oc_signs,confsizes[7])) { caseFound=true; confcase=7; }
  if(!caseFound&&checkCases(confcase8,oc_signs,confsizes[8])) { caseFound=true; confcase=8; }
  if(!caseFound&&checkCases(confcase9,oc_signs,confsizes[9])) { caseFound=true; confcase=9; }
  if(!caseFound&&checkCases(confcase10,oc_signs,confsizes[10])) { caseFound=true; confcase=10; }
  if(!caseFound&&checkCases(confcase11,oc_signs,confsizes[11])) { caseFound=true; confcase=11; }
  if(!caseFound&&checkCases(confcase12,oc_signs,confsizes[12])) { caseFound=true; confcase=12; }
  if(!caseFound&&checkCases(confcase13,oc_signs,confsizes[13])) { caseFound=true; confcase=13; }
  if(!caseFound&&checkCases(confcase14,oc_signs,confsizes[14])) { caseFound=true; confcase=14; }

  if(confcase==6)
  { int a=0; }

  mccase.mcCase=confcase;

  return mccase;
}
//*/

/*
MCCase MCTester::identifyMarchingCubesCase(Cell &values,float iso_val)
{
  float *vals=values.vertexValues;
  char oc_signs[8];
  int confcase=0;
  bool caseFound=false;

  MCCase mccase;

  for(int i=0;i<8;i++) if(vals[i]>iso_val) oc_signs[i] = 0; else oc_signs[i] = 1;
  

  if(!caseFound&&checkCases(confcase0,oc_signs,confsizes[0])) { caseFound=true; confcase=0; }
  if(!caseFound&&checkCases(confcase1,oc_signs,confsizes[1])) { caseFound=true; confcase=1; }
  if(!caseFound&&checkCases(confcase2,oc_signs,confsizes[2])) { caseFound=true; confcase=2; }
  if(!caseFound&&checkCases(confcase3,oc_signs,confsizes[3])) { caseFound=true; confcase=3; }
  if(!caseFound&&checkCases(confcase4,oc_signs,confsizes[4])) { caseFound=true; confcase=4; }
  if(!caseFound&&checkCases(confcase5,oc_signs,confsizes[5])) { caseFound=true; confcase=5; }
  if(!caseFound&&checkCases(confcase6,oc_signs,confsizes[6])) { caseFound=true; confcase=6; }
  if(!caseFound&&checkCases(confcase7,oc_signs,confsizes[7])) { caseFound=true; confcase=7; }
  if(!caseFound&&checkCases(confcase8,oc_signs,confsizes[8])) { caseFound=true; confcase=8; }
  if(!caseFound&&checkCases(confcase9,oc_signs,confsizes[9])) { caseFound=true; confcase=9; }
  if(!caseFound&&checkCases(confcase10,oc_signs,confsizes[10])) { caseFound=true; confcase=10; }
  if(!caseFound&&checkCases(confcase11,oc_signs,confsizes[11])) { caseFound=true; confcase=11; }
  if(!caseFound&&checkCases(confcase12,oc_signs,confsizes[12])) { caseFound=true; confcase=12; }
  if(!caseFound&&checkCases(confcase13,oc_signs,confsizes[13])) { caseFound=true; confcase=13; }
  if(!caseFound&&checkCases(confcase14,oc_signs,confsizes[14])) { caseFound=true; confcase=14; }

  mccase.mcCase=confcase;
  mccase.faceIndex=0;
  mccase.bodyIndex=0;
  getComponentEdges(confcase,oc_signs,mccase.componentEdges);

  if(confcase==0||confcase==1||confcase==2|confcase==5||confcase==8||confcase==9||confcase==11||confcase==14)
    return mccase;

  if(confcase==3)
  {
    int face_idx;
    for(int i=0;i<6;i++)
    {
        if(oc_signs[faces[i][0]]==oc_signs[faces[i][2]]&&
        oc_signs[faces[i][1]]==oc_signs[faces[i][3]]&&
        oc_signs[faces[i][0]]!=oc_signs[faces[i][1]])
        { face_idx=i; break; }
     }

    if(checkFaceSaddle(vals,face_idx,iso_val)==1)
        mccase.faceIndex=2;
    else
        mccase.faceIndex=1;

    return mccase;
  }

  if(confcase==4)
  {
    vector<char> bs_sign;
    computeBodySaddle(vals,iso_val,bs_sign);
    assert(bs_sign.size()>0);
    if(bs_sign[0]==0) mccase.bodyIndex=2;
    if(bs_sign[0]==1) mccase.bodyIndex=1;
    return mccase;
  }

  if(confcase==6)
  {
    int face_idx;
    for(int i=0;i<6;i++)
    {
      if(oc_signs[faces[i][0]]==oc_signs[faces[i][2]]&&
         oc_signs[faces[i][1]]==oc_signs[faces[i][3]]&&
         oc_signs[faces[i][0]]!=oc_signs[faces[i][1]])
      { face_idx=i; break; }
    }

    if(checkFaceSaddle(vals,face_idx,iso_val)==1)
      mccase.faceIndex=2;
    else
    {
      // Check body saddle
      mccase.faceIndex=1;
      
      vector<char> bs_sign;
      computeBodySaddle(vals,iso_val,bs_sign);
      assert(bs_sign.size()>0);
      if(bs_sign[0]==0) mccase.bodyIndex=2;
      if(bs_sign[0]==1) mccase.bodyIndex=1;
    }
    return mccase;
  }

  if(confcase==7)
  {
    int face_idx[3],idx_counter=0;
    int postive_faces=0;
    for(int i=0;i<6;i++)
    {
      if(oc_signs[faces[i][0]]==oc_signs[faces[i][2]]&&
         oc_signs[faces[i][1]]==oc_signs[faces[i][3]]&&
         oc_signs[faces[i][0]]!=oc_signs[faces[i][1]])
        face_idx[idx_counter++]=i;

      if(idx_counter==3)
        break;
    }

    for(int i=0;i<3;i++)
      if(checkFaceSaddle(vals,face_idx[i],iso_val)==1)
        positive_faces++;

    switch(positive_faces)
    {
      case 0: mccases.faceIndex=1;
              break;
      case 1: mccases.faceIndex=2;
              break;
      case 2: mccases.faceIndex=3;
              break;
      case 3:
            {
              mccases.faceIndex=4;
              
              vector<char> bs_sign;
              computeBodySaddle(vals,iso_val,bs_sign);
              assert(bs_sign.size()>0);
              if(bs_sign[0]==0) mccase.bodyIndex=1;
              if(bs_sign[0]==1) mccase.bodyIndex=2;
            }
      break;
    }
    return mccase;
  }

  if(confcase==10||confcase==12)
  {
    int face_idx[2],face_saddle_signs[2],idx_counter=0;

    for(int i=0;i<6;i++)
    {
      if(oc_signs[faces[i][0]]==oc_signs[faces[i][2]]&&
         oc_signs[faces[i][1]]==oc_signs[faces[i][3]]&&
         oc_signs[faces[i][0]]!=oc_signs[faces[i][1]])
        face_idx[idx_counter++]=i;

      if(idx_counter==2)
        break;
    }

    for(int i=0;i<2;i++)
      face_saddle_signs[i]=checkFaceSaddle(vals,face_idx[i],iso_val);

    if(face_saddle_signs[0]==1&&face_saddle_signs[1]==1)
    {
      mccase.faceIndex=1; mccase.bodyIndex=1;
    }
    else if(face_saddle_signs[0]==1&&face_saddle_signs[1]==-1||
            face_saddle_signs[0]==-1&&face_saddle_signs[1]==1)
      mccase.faceIndex=2;
    else
    {
      mccase.faceIndex=1;

      vector<char> bs_sign;
      computeBodySaddle(vals,iso_val,bs_sign);
      assert(bs_sign.size()>0);
      if(bs_sign[0]==0) mccase.bodyIndex=2;
      if(bs_sign[0]==1) mccase.bodyIndex=1;
    }

    return mccase;
  }

  
}
*/

vector<vector<unsigned int> > MCTester::getComponentEdges(int *signs)
{
  int compStart=signs[18],compEnd=signs[19];
  vector<vector<unsigned int> > componentEdges;

  if(compStart!=-1)
  {
    for(int comp=compStart;comp<compEnd;comp++)
    {
      int compSize=components[comp][0];
      vector<unsigned int> component;
      for(int j=1;j<compSize+1;j++)
        component.push_back(components[comp][j]);
      componentEdges.push_back(component);
    }
  }

  return componentEdges;
}