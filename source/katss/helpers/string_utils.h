#ifndef STRING_UTILS_H
#define STRING_UTILS_H


/**
 *  @brief Returns a substring of sequence that is length-characters long starting at 
 *  character start.
 *  
 *  The start parameter defines the starting index of the sequence string
 * 
 *  @note You have to free the substring. Since memory is allocated to store the substring, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 * 
 *  @param sequence The sequence to get the substring from
 *  @param start    The starting index of sequence for substring
 *  @param length   The number of characters in the substring
 * 
 *  @return Pointer to the substring
*/
char *substr(const char *sequence, const int start, const int length);


/**
 *  @brief Searches for the first occurrence of the character `match` in the string pointed
 *  by the argument `str` of size `strsize`.
 * 
 *  This function is useful when searching for a character in a non-null terminating string.
 * 
 *  @param str The string the search in
 *  @param strsize The number of characters in the string
 *  @param match the character to search for
 * 
 *  @return Pointer to the first occurrence of the character `match` in the string `str`, or `NULL`
 *  if the character is not found
*/
char *strnchr(const char *str, size_t strsize, const char match);


/**
 *  @brief Get the basename prefix of a file path.
 * 
 *  This function removes all characters before the last occurrence of the character `/`, and all
 *  characters following the first `'.'`. Assuming `full_path` were to be: `/usr/bin/file.txt`, 
 *  this function would return the string `file`. If the string were to contain no `/` or `.`
 *  characters, then it returns a duplicate of the string.
 * 
 *  @note You have to free the returned string. Since memory is allocated to store the string, 
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 *  @param full_path    The character pointer containing the file path.
 * 
 *  @return char pointer to the basename prefix
*/
char *basename_prefix(const char *full_path);


/**
 *  @brief Concatenate two strings.
 * 
 *  The original strings passed in the argument will remain unaffected, since a new string
 *  containing the combined contents will be returned. 
 * 
 *  @note You have to free the returned string. Since memory is allocated to store the new string,
 *  it is then your responsibility to free the memory when it is no longer in use.
 * 
 *  @param s1   String to concatenate to
 *  @param s2   Appends string to `s1`
 * 
 *  @return char pointer to concatenated string.
*/
char *concat(const char *s1, const char *s2);


/**
 *  @brief Appends the contents of the second string to the end of the first string.
 * 
 *  The function allocates memory for the combined result and updates the first string accordingly.
 *  If the first string is null or empty, it essentially duplicates the contents of `s2` into `s1`.
 *  Similarly, if `s2` is null or empty, then the contents of `s1` will remain the same.
 *
 *  @param s1 The pointer to the first string (modifiable).
 *  @param s2 The second string to append.
 */
void append(char **s1, const char *s2);


/**
 *  @brief Finds the starting index of the first occurrence of a substring in a given string.
 *
 *  This function searches for the first occurrence of the substring specified by `s2` within
 *  the string `s1` and returns the starting index of that occurrence. If the substring is not
 *  found, the function returns -1.
 *
 *  @param s1 The null-terminated string in which the substring is searched.
 *  @param s2 The null-terminated substring to be located within the string `s1`.
 *  @return The starting index of the first occurrence of the substring, or -1 if not found.
 */
int subindx(const char *s1, const char *s2);


/**
 *  @brief Replace every instance of `s2` in `s1` with `X`'s.
 * 
 *  This function searches for every occurrence of the substring `s2` within the string `s1` and
 *  replaces every character that matches `s2` with `X`'s. For example, if `s1` were to be
 *  `Hello wonderful world!` and `s2` were to be `wo`, then `s1` would be modified to look as such:
 *  `Hello XXnderful XXrld!`. If `s2` is not found within `s1`, then `s1` will remain unmodified.
 * 
 *  @param s1 The null-terminated string which will be crossed out
 *  @param s2 The null-termianted substring to be searched for within the string `s1`
*/
void cross_out(char *s1, const char *s2);


/**
 *  @brief Clean a sequence string.
 * 
 *  This function removes trailing newline character returned from rnaf_get(), capitalizes every
 *  nucleotide, and substitutes 'T' and 't' characters with 'U' if specified.
 * 
 *  @param sequence The null-terminated string to clean
 *  @param do_substitute Substitute 'T' and 't' characters with 'U'
*/
void clean_seq(char *sequence, int do_substitute);


/**
 *  @brief Capitalizes every lower case letter in string.
 * 
 *  @param str  String to capitalize
*/
void str_to_upper(char *str);


/**
 *  @brief Converts from DNA alphabet to RNA.
 * 
 *  The function will substitute all 'T' and 't' characters with 'U' and 'u' characters,
 *  respectively.
 * 
 *  @param sequence The sequence to substitute characters
*/
void seq_to_RNA(char *sequence);


void remove_escapes(char *sequence);


#endif
