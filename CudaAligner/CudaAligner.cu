#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <CudaAligner/CudaAligner.h>
#include <CudaAligner/Util.h>

#include <Base/Core.h>

#include <iostream>

struct LocalRowData
{
    int y;
    int row;
    int valBestInCol;
};

__global__ void kernel(int iBlock0, int iBlockRun, int blockRunProgress, int ctColPerThread, int ctRow, int ctCol, int defVal,
    int *row, int *col, int *valBestInRow, int *valBestInCol, char *x, char *y, EndPoint *bests, const Scoring sc)
{
    extern __shared__ LocalRowData _lrd[];

    const int iBlock = iBlock0 + blockIdx.x;

    const int iThread0 = iBlock * blockDim.x;

    LocalRowData *lrd = _lrd - ctColPerThread * iThread0;

    const int iThreadInGrid = iThread0 + threadIdx.x;
    const int iCol0T = ctColPerThread * iThreadInGrid;
    const int iCol1T = ctColPerThread * (iThreadInGrid + 1);

    EndPoint best;
    best.val = INT_MIN;
    best.p = Vec2i{ -1, -1 };

    const int blockOffset = blockRunProgress * (iBlockRun - iBlock) - iBlock - iBlock * 2 * ((int)blockDim.x - 1);

    const int threadCurRowOffset = -2 * (int)threadIdx.x;

    const int iRow0T = blockOffset + threadCurRowOffset;
    const int iRow1T = blockOffset + blockRunProgress + threadCurRowOffset;

    const int iRow1S = min(iRow1T, ctRow);
    
    __syncthreads();

    if (iRow0T < ctRow) {
        for (int iCol = iCol0T; iCol < iCol1T; ++iCol) {
            lrd[iCol].row = row[iCol];
            lrd[iCol].valBestInCol = valBestInCol[iCol];
            lrd[iCol].y = y[iCol];
        }
    }

    __syncthreads();

    for (int iRow = iRow0T; iRow != iRow1T; ++iRow) {
        if (0 <= iRow && iRow < ctRow) {
            int prevColRow = col[iRow - 1];
            // TODO: fix this
            int prevCol = iRow == ctRow - 1 && iThreadInGrid != 0 ? row[iCol0T - 1] : col[iRow];
            int valBestInRowLocal = valBestInRow[iRow];
            char xLocal = x[iRow];
            for (int iCol = iCol0T; iCol < iCol1T; ++iCol) {
                int r = defVal;

                const int prevRow = lrd[iCol].row;
                valBestInRowLocal = max(valBestInRowLocal - sc.k, prevCol - sc.b);
                lrd[iCol].valBestInCol = max(lrd[iCol].valBestInCol - sc.k, prevRow - sc.b);

                r = max(r, valBestInRowLocal);
                r = max(r, lrd[iCol].valBestInCol);

                //r = max(r, prevColRow + sc.match(ca->dX[iRow + rev], ca->dY[iCol + rev]));
                r = max(r, prevColRow + (xLocal == lrd[iCol].y ? sc.mp : sc.mn));

                prevColRow = prevRow;
                lrd[iCol].row = r;
                prevCol = r;

                if (r > best.val && iCol < ctCol) {
                    best.val = r;
                    best.p = Vec2i{ iRow, iCol };
                }
            }
            col[iRow - 1] = prevColRow;
            valBestInRow[iRow] = valBestInRowLocal;
            if (iRow == iRow1S - 1) {
                for (int iCol = iCol0T; iCol < iCol1T; ++iCol) {
                    row[iCol] = lrd[iCol].row;
                    valBestInCol[iCol] = lrd[iCol].valBestInCol;
                }
            }
        }
        __syncthreads();
    }

    // TODO: beter init
    if (iBlockRun == 0 || best.val > bests[iThreadInGrid].val) {
        bests[iThreadInGrid] = best;
    }
} 

int callKernel(int ctRow, int ctCol, char *x, char *y, int defVal, CudaAligner *ca, const Scoring &sc)
{
    const int ctColPerThread = 4;
    const int ctThreadPerFullBlock = 1024;
    const int blockRunProgress = 1024;

    const int ctTotalThreads = (ctCol + ctColPerThread - 1) / ctColPerThread;

    const int ctBlocks = (ctTotalThreads + ctThreadPerFullBlock - 1) / ctThreadPerFullBlock;

    const int ctThreadsPerBlock = (ctTotalThreads + ctBlocks - 1) / ctBlocks;

    const int ctColR = ctBlocks * ctThreadsPerBlock * ctColPerThread;

    const int iLastThreadInBlock = ctThreadsPerBlock - 1;
    const int iLastBlock = ctBlocks - 1;

    const int lastThreadOffset = 2 * iLastThreadInBlock;

    const int ctBlockRuns = (ctRow + blockRunProgress * iLastBlock + iLastBlock + iLastBlock * 2 * iLastThreadInBlock + 2 * iLastThreadInBlock + blockRunProgress - 1) / blockRunProgress;

    for (int iBlockRun = 0; iBlockRun < ctBlockRuns; ++iBlockRun) {
        //const int iBlock0 = min(max(0, (blockRunProgress * iBlockRun - ctRow) / (blockRunProgress + 1)), ctBlocks - 1);
        const int iBlock0 = 0;
        //const int iBlock1 = min(max(0, (blockRunProgress * (iBlockRun + 2)) / (blockRunProgress + 1)), ctBlocks - 1);
        const int iBlock1 = ctBlocks - 1;

        const int ctBlocksTrim = iBlock1 - iBlock0 + 1;

        kernel<<<ctBlocksTrim, ctThreadsPerBlock, ctThreadPerFullBlock * ctColPerThread  * sizeof(LocalRowData)>>>(
            iBlock0, iBlockRun, blockRunProgress, ctColPerThread, ctRow, ctCol, defVal,
            ca->dRow, ca->dCol, ca->dValBestInRow, ca->dValBestInCol, x, y, ca->dBest, sc);

        CUDA_CHECK(cudaGetLastError());

        CUDA_CHECK(cudaDeviceSynchronize());
    }

    return ctTotalThreads;
}
