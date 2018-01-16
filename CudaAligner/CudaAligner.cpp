#include <CudaAligner/CudaAligner.h>

#include <cuda_runtime.h>

#include <iostream>
#include <ctime>

static const int ctThreadsPerBlock = 1024;

extern int callKernel(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int defVal, CudaAligner *ca, const Scoring &sc);

#define CUDA_CHECK(x) cudaCheck(x, __FILE__, __LINE__, #x)

void cudaCheck(cudaError_t status, const char *file, int line, const char *lineStr) {
    if (status != cudaSuccess) {
        cerr << "Cuda error: " << (int)status << " - " << cudaGetErrorString(status) << endl;
        cerr << "  In " << file << " line " << line << endl;
        cerr << "  " << lineStr << endl;
    }
}

EndPoint CudaAligner::_alignCuda(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const vector<char> &x, const vector<char> &y, const Scoring &sc)
{
    const int dir = rev ? -1 : 1;
    const int iRowL = rev ? iRow1 : iRow0 - dir;
    const int iRowH = rev ? iRow0 - dir : iRow1;
    const int iRowL0 = rev ? iRow1 : iRow0;

    const int iColL = rev ? iCol1 : iCol0 - dir;
    const int iColH = rev ? iCol0 - dir : iCol1;
    const int iColL0 = rev ? iCol1 : iCol0;

    col[iRow0 - dir] = row[iCol0 - dir];
    for (int iRow = iRow0; iRow != iRow1 + dir; iRow += dir) {
        col[iRow] = rowStart;
    }

    CUDA_CHECK(cudaMemcpy(dRow + iColL0, row + iColL0, (iColH - iColL0 + 1) * sizeof(*row), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dCol + iRowL, col + iRowL, (iRowH - iRowL + 1) * sizeof(*col), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dValBestInCol + iColL0, valBestInCol.data() + iColL0, (iColH - iColL0 + 1) * sizeof(*dValBestInCol), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dValBestInRow + iRowL0, valBestInRow.data() + iRowL0, (iRowH - iRowL0 + 1) * sizeof(*dValBestInRow), cudaMemcpyHostToDevice));

    //cout << "kernel " << iCol0 << ' ' << iCol1 << endl;
    static int totalTime = 0;
    const int start = clock();

    int ctBest = callKernel(rev, iRow0, iRow1, iCol0, iCol1, rowStart, this, sc);

    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaDeviceSynchronize());

    const int tm = clock() - start;
    totalTime += tm;
    //cout << "time: " << (float)tm / CLOCKS_PER_SEC << ' ' << (float)totalTime / CLOCKS_PER_SEC << endl;

    CUDA_CHECK(cudaMemcpy(row + iColL0, dRow + iColL0, (iColH - iColL0 + 1) * sizeof(*row), cudaMemcpyDeviceToHost));
    //CUDA_CHECK(cudaMemcpy(col + iRowL, dCol + iRowL, (iRowH - iRowL + 1) * sizeof(*col), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(valBestInCol.data() + iColL0, dValBestInCol + iColL0, (iColH - iColL0 + 1) * sizeof(*dValBestInCol), cudaMemcpyDeviceToHost));
    //CUDA_CHECK(cudaMemcpy(valBestInRow.data() + iRowL0, dValBestInRow + iRowL0, (iRowH - iRowL0 + 1) * sizeof(*dValBestInRow), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(bests.data(), dBest, ctBest * sizeof(*dBest), cudaMemcpyDeviceToHost));
    EndPoint best = bests[0];
    for (int i = 1; i < ctBest; ++i) {
        const EndPoint &b = bests[i];
        if (b.val > best.val || (b.val == best.val &&
            (dir * b.p.x < dir * best.p.x || b.p.x == best.p.x && dir * b.p.y < dir * best.p.y))) {
            best = b;
        }
        //best.Add(b);
    }

    return best;
}


void CudaAligner::init(int ctRow, const vector<char> &x, int ctCol, const vector<char> &y)
{
    _col.resize(ctRow + 3);
    col = _col.data() + 1;
    bests.resize(ctThreadsPerBlock);

    // TODO: check sentinels?
    // TODO: allocations alignment?
    cudaMalloc(&dRow, (ctCol + 2) * sizeof(*dRow));
    cudaMalloc(&dCol, (ctRow + 2) * sizeof(*dCol));
    cudaMalloc(&dX, (ctRow + 2) * sizeof(*dX));
    cudaMalloc(&dY, (ctCol + 2) * sizeof(*dY));
    cudaMalloc(&dValBestInRow, (ctRow + 2) * sizeof(*dValBestInRow));
    cudaMalloc(&dValBestInCol, (ctCol + 2) * sizeof(*dValBestInCol));
    // TODO: multiple blocks
    cudaMalloc(&dBest, ctThreadsPerBlock * sizeof(*dBest));

    cudaMemcpy(dX, x.data(), (ctRow + 2) * sizeof(*dX), cudaMemcpyHostToDevice);
    cudaMemcpy(dY, y.data(), (ctCol + 2) * sizeof(*dY), cudaMemcpyHostToDevice);
}

void CudaAligner::destroy()
{
    cudaFree(dBest);
    cudaFree(dValBestInCol);
    cudaFree(dValBestInRow);
    cudaFree(dY);
    cudaFree(dX);
    cudaFree(dCol);
    cudaFree(dRow);
}
