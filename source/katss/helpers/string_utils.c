#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

#include "memory_utils.h"
#include "string_utils.h"

char* substr(const char *sequence, const int start, const int length) {
    char*   substring = (char*) s_malloc((length+1)*sizeof(char));
    int     seq_length = strlen(sequence);
    int     substr_max_length = seq_length - start;

    if(start < 0 || start > seq_length) {
        warning_message("'start' value of %d not valid for sequence '%s'.",start,sequence);
        free(substring);
        return NULL;
    }

    if(!sequence) {
        error_message("Unable to read string %s",sequence);
        return NULL;
    }
    
    if(length > substr_max_length) {
        substring = s_realloc(substring, (substr_max_length+1)*sizeof(char));
        memcpy(substring, &sequence[start], substr_max_length);
        substring[substr_max_length] = '\0';
    } else {
        memcpy(substring, &sequence[start], length);
        substring[length] = '\0';
    }

    return substring;
}


char *strnchr(const char *str, size_t strsize, const char match) {
	for(size_t i=0; i<strsize; i++) {
		if(*str == match) {
			return (char *)str;
		}
		++str;
	}
	return NULL;
}


char *basename_prefix(const char *file_path) {
    char *basename = strrchr(file_path, '/');
    basename++; // Move pointer to remove trailing '/'

    char *ptr = basename;
    int len = 0;
    while(ptr[0] != '.') {
        len++;
        ptr++;
    }

    return substr(basename, 0, len);
}


char *concat(const char *s1, const char *s2) {
    const size_t len1 = strlen(s1);
    const size_t len2 = strlen(s2);

    char *result = s_malloc(len1 + len2 + 1); // +1 for the null-terminator
    if(!result) {
        error_message("Failed to allocate memory for concatenation of '%s' and '%s'.",s1,s2);
    }

    memcpy(result, s1, len1);
    memcpy(result + len1, s2, len2 + 1); // +1 to copy the null-terminator
    return result;
}


void append(char **s1, const char *s2) {
    const size_t len1 = *s1 ? strlen(*s1) : 0;
    const size_t len2 =  s2 ? strlen(s2)  : 0;

    if(len2 == 0) {
        return;     // s2 is NULL or empty, so dont modify contents of s1
    }

    *s1 = s_realloc(*s1,len1 + len2 + 1);

    if(len1 == 0) {
        strcpy(*s1, s2);    // s1 is NULL or empty, so copy contents of s2
    } else {
        strcat(*s1, s2);    // append contents of s2 onto s1
    }
}


int subindx(const char *s1, const char *s2) {
    char *s = strstr(s1, s2);
    if(s) {
        return s - s1;
    } else {
        return -1;
    }
}


void cross_out(char *s1, const char *s2) {
    size_t len_s2 = strlen(s2);
    int indx;
    while( (indx = subindx(s1, s2)) != -1) {
        for(size_t i=indx; i<indx+len_s2; i++) {
            s1[i]='X';
        }
    }
}


void clean_seq(char *sequence, int do_substitute) {
	size_t ln = strlen(sequence)-1;

	if(sequence[ln] == '\n') {  // remove trailing new line character
		sequence[ln] = '\0';
	}

	for(unsigned int i = 0; sequence[i]; i++) {
		sequence[i] = toupper(sequence[i]);
		if(do_substitute && (sequence[i] == 'T' || sequence[i] == 't')) {
			sequence[i] = 'U';
		}
	}
}


void str_to_upper(char *str) {
    if(!str) {
        error_message("Unable to read string %s",str);
    }

    for(int i=0; str[i]; i++) {
        str[i] = toupper(str[i]);
    }
}


void seq_to_RNA(char *sequence) {
    unsigned int i;

    if(!sequence) {
        error_message("Unable to read string %s",sequence);
    }

    for (i = 0; sequence[i]; i++) {
        if (sequence[i] == 'T') {
            sequence[i] = 'U';
        } else if (sequence[i] == 't') {
            sequence[i] = 'u';
        }
    }
}


void remove_escapes(char *str) {

    if(!str) {  // str is NULL
        return;
    }

    size_t ln = strlen(str)-1;

    if(str[ln] == '\n') {  // remove trailing new line character
        str[ln] = '\0';
    }
    
    while(isspace(*str)) {  // move pointer past white space
        ++str;
    }
}
