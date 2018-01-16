#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <CudaAligner/CudaAligner.h>
#include <CudaAligner/Util.h>

#include <Base/Core.h>

struct LocalRowData
{
    int y;
    int row;
    int valBestInCol;
};

__global__ void kernel(int ctColPerThread, int iThread0, bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int defVal,
    int *row, int *col, int *valBestInRow, int *valBestInCol, char *x, char *y, EndPoint *bests, const Scoring sc)
{
    extern __shared__ LocalRowData _lrd[];

    const int dir = rev ? -1 : 1;

    LocalRowData *lrd = _lrd - (!rev ? iCol0 + dir * ctColPerThread * iThread0 : iCol0 + dir * ctColPerThread * (iThread0 + blockDim.x) - dir);

    iRow1 += dir;
    iCol1 += dir;

    //const int iCol0T = iCol0 + (iCol1 - iCol0) * threadIdx.x / blockDim.x;
    //const int iCol1T = iCol0 + (iCol1 - iCol0) * (threadIdx.x + 1) / blockDim.x;
    //printf("%d %d %d %d %d %d\n", iCol0T, iCol1T, iCol1, iCol1 - iCol0, (iCol1 - iCol0) * (threadIdx.x + 1), (iCol1 - iCol0) * (threadIdx.x + 1) / blockDim.x);

    const int iCol0T = iCol0 + dir * ctColPerThread * (iThread0 + threadIdx.x);
    int iCol1T = iCol0 + dir * ctColPerThread * (iThread0 + threadIdx.x + 1);
    iCol1T = !rev ? min(iCol1T, iCol1) : max(iCol1T, iCol1);

    EndPoint best;
    best.val = INT_MIN;
    best.p = Vec2i{ -1, -1 };

    int iRun = 0;
    const int threadRowDiff = 2;
    const int threadRowProgress = threadRowDiff - 1;
    const int threadCurRowOffset = dir * (threadRowDiff * (iRun - threadIdx.x) - threadIdx.x);
    const int lastThreadCurRowOffset = dir * (threadRowDiff * (iRun - (blockDim.x - 1)) - (blockDim.x - 1));

    const int iRow0T = iRow0 + threadCurRowOffset;
    const int iRow1T = iRow1 + threadCurRowOffset - lastThreadCurRowOffset;

    __syncthreads();

    //printf("%d %d\n", iCol0T, iCol1T);

    //for (int iCol = iCol0T; iCol != iCol1T; iCol += dir) {
    //    lrd[iCol].row = row[iCol];
    //    lrd[iCol].valBestInCol = valBestInCol[iCol];
    //    lrd[iCol].y = y[iCol + rev];
    //}

    //__syncthreads();

    //const bool cnd = iRow0 == 8 && iRow1 == 9;

    for (int iRow = iRow0T; iRow != iRow1T; iRow += dir) {
        if (!rev && iRow0 <= iRow && iRow < iRow1 || rev && iRow0 >= iRow && iRow > iRow1) {
            int prevColRow = col[iRow - dir];
            // TODO: fix this
            int prevCol = iRow == iRow1 - dir && threadIdx.x + iThread0 != 0 ? row[iCol0T - dir] : col[iRow];
            int valBestInRowLocal = valBestInRow[iRow];
            char xLocal = x[iRow + rev];
            for (int iCol = iCol0T; iCol != iCol1T; iCol += dir) {
                int r = defVal;

                const int prevRow = row[iCol];
                valBestInRowLocal = max(valBestInRowLocal - sc.k, prevCol - sc.b);
                valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, prevRow - sc.b);

                r = max(r, valBestInRowLocal);
                r = max(r, valBestInCol[iCol]);

                //r = max(r, prevColRow + sc.match(ca->dX[iRow + rev], ca->dY[iCol + rev]));
                r = max(r, prevColRow + (xLocal == y[iCol + rev] ? sc.mp : sc.mn));

                //if (cnd) {
                //    printf("%d %d %d %d %d\n", iCol, r, prevRow, prevColRow, prevCol);
                //}

                prevColRow = prevRow;
                row[iCol] = r;
                prevCol = r;

                if (r > best.val) {
                    best.val = r;
                    best.p = Vec2i{ iRow, iCol };
                }
            }
            col[iRow - dir] = prevColRow;
            valBestInRow[iRow] = valBestInRowLocal;
        }
        ++iRun;
        __syncthreads();
    }

    bests[iThread0 + threadIdx.x] = best;
}

int callKernel(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int defVal, CudaAligner *ca, const Scoring &sc)
{
    //const int ctCol = abs(iCol1 - iCol0) + 1;
    //int ctThread = min(1024, ctCol);
    //const int ctColPerThread = ctCol / ctThread;
    //const int ctColExtra = ctCol - ctColPerThread * ctThread;
    //kernel<<<1, ctThread>>>(ctColPerThread, ctColExtra, rev, iRow0, iRow1, iCol0, iCol1, defVal,
    //    ca->dRow, ca->dCol, ca->dValBestInRow, ca->dValBestInCol, ca->dX, ca->dY, ca->dBest, sc);
    //return ctThread;

    const int ctColPerThread = 4;
    const int ctThreadPerFullBlock = 1024;

    const int ctCol = abs(iCol1 - iCol0) + 1;

    //const int ctColPerFullBlock = ctThreadPerFullBlock * ctColPerThread;

    const int ctTotalThreads = (ctCol + ctColPerThread - 1) / ctColPerThread;

    const int ctBlocks = (ctTotalThreads + ctThreadPerFullBlock - 1) / ctThreadPerFullBlock;

    for (int iBlock = 0; iBlock < ctBlocks; ++iBlock) {
        const int iThread0 = ctTotalThreads * iBlock / ctBlocks;
        const int iThread1 = ctTotalThreads * (iBlock + 1) / ctBlocks;

        const int ctThread = iThread1 - iThread0;

        kernel<<<1, ctThread, ctThreadPerFullBlock  * ctColPerThread  * sizeof(LocalRowData)>>>(
            ctColPerThread, iThread0, rev, iRow0, iRow1, iCol0, iCol1, defVal,
            ca->dRow, ca->dCol, ca->dValBestInRow, ca->dValBestInCol, ca->dX, ca->dY, ca->dBest, sc);

        CUDA_CHECK(cudaGetLastError());

        CUDA_CHECK(cudaDeviceSynchronize());
    }

    return ctTotalThreads;
}
