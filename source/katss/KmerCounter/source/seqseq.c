/* 
 * seqseq.c - Locate the first occurrence of a nucleotide sequence in another sequence.
 *
 * The 'seqseq' function is based on the strstr implementation from the GNU C Library (glibc).
 *
 * Copyright (C) 1994-2024 Free Software Foundation, Inc.
 * This file is part of the GNU C Library.
 *
 * The GNU C Library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * The GNU C Library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with the GNU C Library; if not, see
 * <https://www.gnu.org/licenses/>.
 *
 * Modifications and additions for nucleotide sequence matching:
 * (C) 2024 Francisco F. Cavazos
 * - Added clean_nt function for case-insensitive and U-to-T conversion
 * - Modified seqseq function to use clean_nt for nucleotide sequences
 * - Added seqncmp and seqchr functions for nucleotide sequence comparison
 *
 * This modified code retains the LGPL license of the original work.
 */

#include <stddef.h>
#include <stdint.h>
#include <string.h>

static inline int clean_nt(const char c1);
static inline int seqcmp(const char *seq, const char *pat);
static inline int seqncmp(const unsigned char *seq, const unsigned char *pat, size_t len);
static inline int seqcmpa(const char *seq, const char *pat);
static inline char *seqchr(const char *seq, int nt);

/*==================================================================================================
|                                    Sequence Search Functions                                     |
==================================================================================================*/
static inline char *
seqseq2(const unsigned char *hs, const unsigned char *ne)
{
	uint32_t h1 = (clean_nt(ne[0]) << 16) | clean_nt(ne[1]);
	uint32_t h2 = 0;
	for(int c = clean_nt(hs[0]); h1 != h2 && c != 0; c = clean_nt(*++hs))
		h2 = (h2 << 16) | c;
	return h1 == h2 ? (char *)hs - 2 : NULL;
}


static inline char *
seqseq3(const unsigned char *hs, const unsigned char *ne)
{
	uint32_t h1 = ((uint32_t)clean_nt(ne[0]) << 24) | (clean_nt(ne[1]) << 16) | (clean_nt(ne[2]) << 8);
	uint32_t h2 = 0;
	for(int c = clean_nt(hs[0]); h1 != h2 && c != 0; c = clean_nt(*++hs))
		h2 = (h2 | c) << 8;
	return h1 == h2 ? (char *)hs - 3 : NULL;
}

/* Hash character pairs so a small shift table can be used.  All bits of
   p[0] are included, but not all bits from p[-1].  So if two equal hashes
   match on p[-1], p[0] matches too.  Hash collisions are harmless and result
   in smaller shifts.  */
#define hash2(p) (((size_t)clean_nt((p)[0]) - ((size_t)clean_nt((p)[-1]) << 3)) % sizeof(shift))


char *
seqseq(const char *seq, const char *pat)
{
	register const unsigned char *hs = (const unsigned char *)seq;
	register const unsigned char *ne = (const unsigned char *)pat;

	/* Handle short needle special cases first.  */
	if(ne[0] == '\0')
		return (char *)hs;
	hs = (const unsigned char *)seqchr((const char*)hs, ne[0]);
	if(hs == NULL || ne[1] == '\0')
		return (char*)hs;
	if(ne[2] == '\0')
		return seqseq2(hs, ne);
	if(ne[3] == '\0')
		return seqseq3(hs, ne);

	/* Ensure haystack length is at least as long as needle length.
	   Since a match may occur early on in a huge haystack, use strnlen
	   and read ahead a few cachelines for improved performance.  */
	size_t ne_len = strlen((const char*)ne);
	size_t hs_len = strnlen((const char*)hs, ne_len | 512);
	if(hs_len < ne_len)
		return NULL;

	/* Check whether we have a match.  This improves performance since we
	   avoid initialization overheads.  */
	if(seqncmp(hs, ne, ne_len))
		return (char *) hs;

	if(ne_len > 256)
		return NULL;

	const unsigned char *end = hs + hs_len - ne_len;
	uint8_t shift[256];
	register size_t tmp, shift1;
	register size_t m1 = ne_len - 1;
	register size_t offset = 0;

	/* Initialize bad character shift hash table.  */
	memset(shift, 0, sizeof (shift));
	for(size_t i = 1; i < m1; i++)
		shift[hash2(ne + i)] = i;
	/* Shift1 is the amount we can skip after matching the hash of the
	   needle end but not the full needle.  */
	shift1 = m1 - shift[hash2(ne + m1)];
	shift[hash2(ne + m1)] = m1;

	while(1) {
		if(hs > end) {
			end += strnlen((const char*)end + m1 + 1, 2048);
			if (hs > end)
				return NULL;
		}

		/* Skip past character pairs not in the needle.  */
		do {
			hs += m1;
			tmp = shift[hash2(hs)];
		} while(tmp == 0 && hs <= end);

		/* If the match is not at the end of the needle, shift to the end
		   and continue until we match the hash of the needle end.  */
		hs -= tmp;
		if(tmp < m1)
			continue;

		/* Hash of the last 2 characters matches.  If the needle is long,
		   try to quickly filter out mismatches.  */
		if(m1 < 15 || seqncmp(hs + offset, ne + offset, 8)) {
			if(seqncmp(hs, ne, m1))
				return (void *)hs;

			/* Adjust filter offset when it doesn't find the mismatch.  */
			offset = (offset >= 8 ? offset : m1) - 8;
		}

		/* Skip based on matching the hash of the needle end.  */
		hs += shift1;
	}
}


char *
seqlseq(const char *seq, const char *pat)
{
	register char *ret = seqseq(seq, pat);
	if(ret == NULL)
		return NULL;

	while(ret > seq && *(ret - 1) != '\n') ret--;
	return ret;
}


char *
seqseqa(const char *seq, const char *pat)
{
	register int c, sc;

	/* Check pat > 0 */
	if((c = clean_nt(*pat++)) == 0)
		return (char *)seq;

	do {
		do {
			switch((sc = clean_nt(*seq++))) {
			case '\0': return NULL;
			case '>': while(*seq++ != '\n' && *seq); break;
			}
		} while(sc != c);
	} while(!seqcmpa(seq, pat));
	seq--;

	return (char *)seq;
}


char *
seqlseqa(const char *seq, const char *pat)
{
	register int c, sc;
	register char *l = NULL;

	/* Check pat > 0 */
	if((c = clean_nt(*pat++)) == 0)
		return (char *)seq;

	do {
		do {
			switch((sc = clean_nt(*seq++))) {
			case '\0':
				return NULL;
			case '>':
				while(*seq++ != '\n' && *seq);
				l = (char *)seq;
				break;
			}
		} while(sc != c);
	} while(!seqcmpa(seq, pat));

	return l;
}


char *
seqseqq(const char *seq, const char *pat)
{
	register int c, sc, cnt;

	/* Check pat > 0 */
	if((c = clean_nt(*pat++)) == 0)
		return (char *)seq;

	do {
		do {
			switch((sc = clean_nt(*seq++))) {
			case '\0':
				return NULL;
			case '@':
				while(*seq++ != '\n' && *seq);
				break;
			case '+': 
				cnt = 0;
				while(cnt < 3 && *seq) if(*seq++ == '\n') cnt++;
				break;
			}
		} while(sc != c);
	} while(!seqcmp(seq, pat));
	seq--;

	return (char *)seq;
}


char *
seqlseqq(const char *seq, const char *pat)
{
	register int c, sc, cnt;
	register char *l;

	/* Check pat > 0 */
	if((c = clean_nt(*pat++)) == 0)
		return (char *)seq;

	do {
		do {
			switch((sc = clean_nt(*seq++))) {
			case '\0':
				return NULL;
			case '@':
				while(*seq++ != '\n' && *seq);
				l = (char *)seq;
				break;
			case '+': 
				cnt = 0;
				while(cnt < 3 && *seq) if(*seq++ == '\n') cnt++;
				l = (char *)seq;
				break;
			}
		} while(sc != c);
	} while(!seqcmp(seq, pat));

	return l;
}

/*==================================================================================================
|                                       Non-search Functions                                       |
==================================================================================================*/
/**
 * @brief Converts all lowercase letters to uppercase & converts 'U' -> 'T'. Meant
 * to clean up nucleotide sequences.
 * 
 * @param c1 Nucleotide character to clean
 * @return int New nucleotide character
 */
static inline int
clean_nt(const char c1)
{
	register int ret = (int)c1;

	/* Check if lowercase letter */
	if(97 <= c1 && c1 <= 122)
		ret = c1 - 32;

	/* Check if U, convert to T */
	if(ret == 85)
		ret--;

	return ret;
}


/**
 * @brief Compare a sequence with a pattern. 
 * 
 * This function determines if `seq` string begins with the `pat` string.
 * Moves the sequence pointer while comparing.
 * 
 * @param seq Sequence that is being searched
 * @param pat Pattern to check if it begin in `seq`
 * @return int `0` if seq does not begin with pat. `1` if `seq` begins with
 * `pat`.
 */
static inline int
seqcmp(const char *seq, const char *pat)
{
	while(*pat)
		if(! *seq || clean_nt(*pat++) != clean_nt(*seq++))
			return 0;
	return 1;
}


/**
 * @brief Compare a sequence with a pattern. 
 * 
 * This function determines if `seq` string begins with the `pat` string.
 * 
 * @param seq Sequence to compare
 * @param pat Pattern to compare to
 * @param len Compare up to `len` characters
 * @return int 
 */
static inline int
seqncmp(const unsigned char *seq, const unsigned char *pat, size_t len)
{
	while(*pat && len--)
		if(! *seq || clean_nt(*pat++) != clean_nt(*seq++))
			return 0;
	return 1;
}


/**
 * @brief Compare a sequence with a pattern, ignoring newline characters in
 * `seq`.
 * 
 * This function determines if `seq` string begins with the `pat` string.
 * Assumes that both `seq` and `pat` are null terminated strings.
 * 
 * @param seq Sequence that is being searched
 * @param pat Pattern to check if it begin in `seq`
 * @return int `0` if seq does not begin with pat. `1` if `seq` begins with
 * `pat`.
 */
static inline int
seqcmpa(const char *seq, const char *pat)
{
	while(*pat) {
		if(*seq == '\0')
			return 0;
		if(*seq == '\n') {
			seq++;
			continue;
		}
		if(clean_nt(*pat++) != clean_nt(*seq++))
			return 0;
	}
	return 1;
}


/**
 * @brief Locate the first occurrence of a nucleotide in a sequence
 * 
 * @param seq Sequence to test on
 * @param nt  nt to search for
 * @return char* Location of character, NULL if nucleotide not in seq.
 */
static inline char *
seqchr(const char *seq, int nt)
{
	register int snt = clean_nt((char)nt);

	while(clean_nt(*seq) != snt)
		if(*seq++ == '\0')
			return NULL;
	return (char *)seq;
}
