/***** shulen.c **************************************************************
 * Description: Compute the lengths of shortest unique substrings.
 * Author: Bernhard Haubold, bernhard.haubold@fh-weihenstephan.de
 * File created on Wed Dec  6 16:12:32 2006.
 * Modified by Mirjana Domazet-Loso, September 22, 2008
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA.
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <divsufsort.h>
#include "commonSC.h"
#include "eprintf.h"
#include "sequenceData.h"
#include "shulen.h"

Int64 *getDivSortSa(Sequence *seq) {
  sauchar_t *t; // The text
  saidx_t *sa;  // The sa
  Int64 i, n, *sa2;

  // Get text
  t = (sauchar_t *)seq->seq;
  // Get length of text
  n = seq->len;
  // Calculate sa. Alfy assumes that the sa starts at sa[1] rather
  // than sa[0]. This is a hold over from using the deep shallow
  // library.
  sa = (saidx_t *)emalloc((size_t)(n+1) * sizeof(saidx_t));
  if (divsufsort(t, sa+1, (saidx_t)n) != 0) {
    printf("ERROR[getDivSortSa]: suffix sorting failed.\n");
    exit(-1);
  }
  sa2 = (Int64 *)malloc((n+1) * sizeof(Int64));
  for (i = 0; i <= n; i++)
    sa2[i] = sa[i];
  free(sa);

  return sa2;
}

/* getLcp calculates the lcp array using Kasais' algorithm. */
Int64 *getLcp(Sequence *seq, Int64 *sa){
  char *t = seq->seq;
  Int64 n = seq->len;
  Int64 i, j, k, l;
  Int64 *isa = (Int64 *)malloc((n+1) * sizeof(Int64));
  Int64 *lcp = (Int64 *)malloc((n+1) * sizeof(Int64));
  
  for (i = 0; i < n+1; i++)
    isa[sa[i]] = i;
  lcp[0] = -1;
  l = 0;
  for (i = 0; i < n+1; i++) {
    j = isa[i];
    if (j == 0)
      continue;
    k = sa[j-1];
    while (t[i+l] == t[k+l])
      l++;
    lcp[j] = l;
    l--;
    if (l < 0)
      l = 0;
  }
  
  return lcp;
}
