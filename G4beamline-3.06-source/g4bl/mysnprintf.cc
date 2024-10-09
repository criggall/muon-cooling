//	mysnprintf.cc -- workaround for missing snprintf

#include <string.h>
#include <stdio.h>
#include "mysnprintf.hh"
#include "globals.hh"
#include "BLAssert.hh"

int myvsnprintf(char *dest, size_t dest_size, const char *format, va_list ap)
{
	char tmp[5000];

	BLAssert(dest_size < 5000);

	int i = vsprintf(tmp,format,ap);
	strncpy(dest,tmp,dest_size);
	if((size_t)i >= dest_size) dest[dest_size-1] = '\0';

	return i;
}

int mysnprintf(char *dest, size_t dest_size, const char *format, ...)
{
        va_list args;
        va_start(args, format);

        int i = myvsnprintf(dest,dest_size,format,args);

        va_end(args);
        return i;
}

