#include <cuda_runtime.h>

#include <iostream>

using namespace std;

void cudaCheck(cudaError_t status, const char *file, int line, const char *lineStr)
{
    if (status != cudaSuccess) {
        cerr << "Cuda error: " << (int)status << " - " << cudaGetErrorString(status) << endl;
        cerr << "  In " << file << " line " << line << endl;
        cerr << "  " << lineStr << endl;
    }
}