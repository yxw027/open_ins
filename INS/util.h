#ifndef VBASE_UTIL_H
#define VBASE_UTIL_H
#include <stdint.h>
#ifndef MAXFIELD
#define MAXFIELD 100
#endif
#ifndef PI
#define PI (3.1415926535897932384626433832795)
#endif


int8_t is_equal_ignore_case(const char* s1, const char* s2);
int8_t parse_fields(char* const buffer, char** val, int* size);



#endif

