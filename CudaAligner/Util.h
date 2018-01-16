#pragma once

#include <cuda_runtime.h>

#define CUDA_CHECK(x) cudaCheck(x, __FILE__, __LINE__, #x)

void cudaCheck(cudaError_t status, const char *file, int line, const char *lineStr);
