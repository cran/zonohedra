#include <windows.h>

static nanotime_t get_nanotime(void) {

    static  LARGE_INTEGER   frequency = {0,0} ;

    if( frequency.QuadPart == 0 )    ::QueryPerformanceFrequency(&frequency);

    LARGE_INTEGER time_var;

    ::QueryPerformanceCounter(&time_var);

    /* Convert to nanoseconds */
    return (nanotime_t)(1.0e9 * time_var.QuadPart / frequency.QuadPart);
}
