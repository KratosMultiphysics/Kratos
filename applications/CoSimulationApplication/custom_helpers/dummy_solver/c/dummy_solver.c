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


void ImportData(const char* pIdentifier)
{
    int incoming_data_size = 1;

    double**
    values_raw = (double**)malloc(sizeof(double*)*1);
    values_raw[0]=(double*)malloc(sizeof(double)*incoming_data_size);

    CoSimIO_ImportData("aaa", pIdentifier, &incoming_data_size, values_raw);

    free(values_raw[1]);
    free(values_raw);
}

void ExportData(const char* pIdentifier)
{

}

void ImportMesh(const char* pIdentifier)
{

}

void ExportMesh(const char* pIdentifier)
{

}


int main()
{
    CoSimIO_Connect("aaa", "ccc");

    CoSimIO_Register("aaa", "ImportData", &ImportData);
    CoSimIO_Register("aaa", "ExportData", &ExportData);
    CoSimIO_Register("aaa", "ImportMesh", &ImportMesh);
    CoSimIO_Register("aaa", "ExportMesh", &ExportMesh);

    CoSimIO_Disconnect("aaa");

    // printf() displays the string inside quotation
    printf("Hello, World!");
    return 0;
}