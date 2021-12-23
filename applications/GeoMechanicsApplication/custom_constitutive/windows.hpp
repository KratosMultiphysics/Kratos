// including windows.h adds some undesirable macros. Here, min and max macros are deactivated.
#ifdef NOMINMAX
#include <windows.h>
#else
#define NOMINMAX
#include <windows.h>
#undef NOMINMAX
#endif