// including Windows.h adds some undesirable macros. Here, min and max macros are deactivated.
#ifdef NOMINMAX
#include <Windows.h>
#else
#define NOMINMAX
#include <Windows.h>
#undef NOMINMAX
#endif