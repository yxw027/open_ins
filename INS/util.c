#include <assert.h>
#include <string.h>
#include "util.h"


int8_t is_equal_ignore_case(const char* s1, const char* s2)
{
	assert(s1 != NULL && s2 != NULL);

	do {
		char c1 = *s1, c2 = *s2;
		if (c1 >= 'A' && c1 <= 'Z') {
			c1 -= ('A' - 'a');
		}
		if (c2 >= 'A' && c2 <= 'Z') {
			c2 -= ('A' - 'a');
		}
		if (c1 != c2) {
			return 0;
		}
		s1++, s2++;
	} while (*s1 != 0 || *s2 != 0);

	return 1;
}

int8_t  parse_fields(char* const buffer, char** val, int* size)
{
	int8_t ret = -1;
	char* p, *q;
	int n = 0;

	/* parse fields */
	for (p = buffer; *p && n < MAXFIELD; p = q + 1) {
		if ((q = strchr(p, ',')) || (q = strchr(p, '*'))) {
			val[n++] = p; 
			*q = '\0';
		}
		else
		{
			val[n++] = p;
			break;
		}
			
	}
	*size = n;
	return 1;

}

