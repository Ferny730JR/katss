/* seqfstrerror.c - seqf functions to get a readable error output
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include <errno.h>
#include <string.h>

#include "seqf_core.h"

static const char seqf_err_msg[7][60] = {
	"No error",
	"Mutex failed to initialize",
	"Invalid mode passed to seqfopen",
	"Read failed, could not determine type of file",
	"Read failed, sequence is larger than input buffer",
	"Out of memory",
	"gets failed, sequence is larger than passed buffer"
};

static const char seqf_undeferr[19] = "Unrecognized error";

int
seqfstrerror_r(int _rnaferrno, char *buffer, size_t bufsize)
{
	if(_rnaferrno == 1)
		return strerror_r(errno, buffer, bufsize);
	if(0 <= _rnaferrno && _rnaferrno <= 7)
		strncpy(buffer, seqf_err_msg[_rnaferrno], bufsize);
	else
		strncpy(buffer, seqf_undeferr, bufsize);
	if(bufsize < strlen(seqf_err_msg[_rnaferrno])) // not enough space
		return 1;
	return 0;
}

const char *
seqfstrerror(int _rnaferrno)
{
	if(_rnaferrno == 1)
		return strerror(errno);
	if(0 <= _rnaferrno && _rnaferrno <= 7)
		return seqf_err_msg[_rnaferrno];
	return seqf_undeferr;
}
