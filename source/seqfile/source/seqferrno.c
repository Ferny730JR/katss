/* seqferrno.c - seqf's errno implementation
 * 
 * Copyright (c) Francisco F. Cavazos 2024
 * Subject to the MIT License
 */

#include "seqf_core.h"

_Thread_local int seqferrno_;

int *
seqfgeterrno(void)
{
	return &seqferrno_;
}
