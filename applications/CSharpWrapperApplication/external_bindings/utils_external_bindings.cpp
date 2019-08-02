#define EXPORT __declspec(dllexport)


extern "C" {
EXPORT void __stdcall Utils_DeleteDoubleArray(double *array) {
    delete[] array;
}
}
