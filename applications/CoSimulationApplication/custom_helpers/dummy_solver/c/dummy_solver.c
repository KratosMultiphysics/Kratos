#include <stdio.h>
#include <stdlib.h>
#include "co_sim_c_io.h"

void abc() {
    printf("Hello, abc!");
}

 int puts(const char *s)
 {
    printf("Called the puts function!");
    return 0;
 }



int main()
{
    CoSimIO_Connect("aaa", "ccc");


    int incoming_data_size = 1;

    double**
    values_raw = (double**)malloc(sizeof(double*)*1);
    values_raw[0]=(double*)malloc(sizeof(double)*incoming_data_size);

    CoSimIO_ImportData("aaa", "bbb", &incoming_data_size, values_raw);

    free(values_raw[1]);
    free(values_raw);

    // printf() displays the string inside quotation
    printf("Hello, World!");
    return 0;
}