#include <stdio.h>
#include "co_sim_c_io.h"

void abc() {
    printf("Hello, abc!");
}

 int puts(const char *s)
 {
    printf("Called the puts function!");
 }



int main()
{
    CoSimIO_Connect("aaa");

    // printf() displays the string inside quotation
    printf("Hello, World!");
    return 0;
}