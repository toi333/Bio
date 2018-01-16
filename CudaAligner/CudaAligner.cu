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

__global__ void kernel(int iBlockRun, int blockRunProgress, int ctColPerThread, bool rev,
    int iRow0, int iRow1, int iCol0, int iCol1, int iCol1Real, int defVal,
    int *row, int *col, int *valBestInRow, int *valBestInCol, char *x, char *y, EndPoint *bests, const Scoring sc)
{
    extern __shared__ LocalRowData _lrd[];

    const int dir = rev ? -1 : 1;

    const int iThread0 = blockIdx.x * blockDim.x;

    LocalRowData *lrd = _lrd - (!rev ? iCol0 + dir * ctColPerThread * iThread0 : iCol0 + dir * ctColPerThread * (iThread0 + (int)blockDim.x) - dir);

    iRow1 += dir;
    iCol1 += dir;

    //const int iCol0T = iCol0 + (iCol1 - iCol0) * threadIdx.x / blockDim.x;
    //const int iCol1T = iCol0 + (iCol1 - iCol0) * (threadIdx.x + 1) / blockDim.x;
    //printf("%d %d %d %d %d %d\n", iCol0T, iCol1T, iCol1, iCol1 - iCol0, (iCol1 - iCol0) * (threadIdx.x + 1), (iCol1 - iCol0) * (threadIdx.x + 1) / blockDim.x);

    const int iCol0T = iCol0 + dir * ctColPerThread * (iThread0 + (int)threadIdx.x);
    const int iCol1T = iCol0 + dir * ctColPerThread * (iThread0 + (int)threadIdx.x + 1);

    EndPoint best;
    best.val = INT_MIN;
    best.p = Vec2i{ -1, -1 };

    const int blockOffset = dir * (blockRunProgress * (iBlockRun - (int)blockIdx.x) - (int)blockIdx.x);

    const int iRow0B = iRow0 + blockOffset;
    const int iRow1B = iRow0B + dir * blockRunProgress;

    const int threadRowProgress = 1;
    const int threadCurRowOffset = dir * (threadRowProgress * (-(int)threadIdx.x) - (int)threadIdx.x);
    const int lastThreadCurRowOffset = dir * (threadRowProgress * (-((int)blockDim.x - 1)) - ((int)blockDim.x - 1));

    const int iRow0T = iRow0B + threadCurRowOffset;
    const int iRow1T = iRow1B + threadCurRowOffset - lastThreadCurRowOffset;

    const int iRow0S = !rev ? max(iRow0, iRow0B) : min(iRow0, iRow0B);
    const int iRow1S = !rev ? min(iRow1, iRow1B) : max(iRow1, iRow1B);

    __syncthreads();

    //printf("%d %d\n", iCol0T, iCol1T);

    if (dir * iRow0T < dir * iRow1) {
        for (int iCol = iCol0T; iCol != iCol1T; iCol += dir) {
            lrd[iCol].row = row[iCol];
            lrd[iCol].valBestInCol = valBestInCol[iCol];
            lrd[iCol].y = y[iCol + rev];
        }
    }

    __syncthreads();

    //const bool cnd = iRow0 == 8 && iRow1 == 9;

    for (int iRow = iRow0T; iRow != iRow1T; iRow += dir) {
        if (dir * iRow0S <= dir * iRow && dir * iRow < dir * iRow1S) {
            int prevColRow = col[iRow - dir];
            // TODO: fix this
            const bool end = iRow == iRow1S - dir || iRow == iRow1S - dir;
            int prevCol = end && threadIdx.x + iThread0 != 0 ? row[iCol0T - dir] : col[iRow];
            int valBestInRowLocal = valBestInRow[iRow];
            char xLocal = x[iRow + rev];
            for (int iCol = iCol0T; iCol != iCol1T; iCol += dir) {
                int r = defVal;

                const int prevRow = lrd[iCol].row;
                valBestInRowLocal = max(valBestInRowLocal - sc.k, prevCol - sc.b);
                lrd[iCol].valBestInCol = max(lrd[iCol].valBestInCol - sc.k, prevRow - sc.b);

                r = max(r, valBestInRowLocal);
                r = max(r, lrd[iCol].valBestInCol);

                //r = max(r, prevColRow + sc.match(ca->dX[iRow + rev], ca->dY[iCol + rev]));
                r = max(r, prevColRow + (xLocal == lrd[iCol].y ? sc.mp : sc.mn));

                //if (cnd) {
                //    printf("%d %d %d %d %d\n", iCol, r, prevRow, prevColRow, prevCol);
                //}

                prevColRow = prevRow;
                lrd[iCol].row = r;
                prevCol = r;

                if (r > best.val && dir * iCol <= dir * iCol1Real) {
                    best.val = r;
                    best.p = Vec2i{ iRow, iCol };
                }
            }
            col[iRow - dir] = prevColRow;
            valBestInRow[iRow] = valBestInRowLocal;
            if (end) {
                for (int iCol = iCol0T; iCol != iCol1T; iCol += dir) {
                    row[iCol] = lrd[iCol].row;
                    valBestInCol[iCol] = lrd[iCol].valBestInCol;
                }
            }
        }
        __syncthreads();
    }

    // TODO: beter init
    if (iBlockRun == 0 || best.val > bests[iThread0 + threadIdx.x].val) {
        bests[iThread0 + threadIdx.x] = best;
    }
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
    const int blockRunProgress = 1024;

    const int ctCol = abs(iCol1 - iCol0) + 1;

    //const int ctColPerFullBlock = ctThreadPerFullBlock * ctColPerThread;

    const int ctTotalThreads = (ctCol + ctColPerThread - 1) / ctColPerThread;

    const int ctBlocks = (ctTotalThreads + ctThreadPerFullBlock - 1) / ctThreadPerFullBlock;

    const int ctThreadsPerBlock = (ctTotalThreads + ctBlocks - 1) / ctBlocks;

    const int dir = rev ? -1 : 1;
    const int iCol1R = iCol0 + dir * (ctBlocks * ctThreadsPerBlock * ctColPerThread - 1);

    const int iLastBlock = ctBlocks - 1;

    const int ctBlockRuns = ((abs(iRow1 - iRow0) + 1 + iLastBlock) + blockRunProgress - 1) / blockRunProgress + iLastBlock;

    for (int iBlockRun = 0; iBlockRun < ctBlockRuns; ++iBlockRun) {
        kernel<<<ctBlocks, ctThreadsPerBlock, ctThreadPerFullBlock  * ctColPerThread  * sizeof(LocalRowData)>>>(
            iBlockRun, blockRunProgress, ctColPerThread, rev, iRow0, iRow1, iCol0, iCol1R, iCol1, defVal,
            ca->dRow, ca->dCol, ca->dValBestInRow, ca->dValBestInCol, ca->dX, ca->dY, ca->dBest, sc);

        CUDA_CHECK(cudaGetLastError());

        CUDA_CHECK(cudaDeviceSynchronize());
    }

    return ctTotalThreads;
}
