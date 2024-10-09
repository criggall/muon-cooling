//	mysnprintf.hh -- workaround for missing snprintf()

#include <stdarg.h>

#define snprintf mysnprintf
#define vsnprintf myvsnprintf

int mysnprintf(char *dest, size_t dest_size, const char *format, ...);
int myvsnprintf(char *dest, size_t dest_size, const char *format, va_list ap);


