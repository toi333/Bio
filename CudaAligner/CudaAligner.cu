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

__global__ void kernel(int iBlockRun, int blockRunProgress, int ctColPerThread, int ctRow, int ctCol, int ctColReal, int defVal,
    int *row, int *col, int *valBestInRow, int *valBestInCol, char *x, char *y, EndPoint *bests, const Scoring sc)
{
    extern __shared__ LocalRowData _lrd[];

    const int iThread0 = blockIdx.x * blockDim.x;

    LocalRowData *lrd = _lrd - ctColPerThread * iThread0;

    const int iThreadInGrid = iThread0 + threadIdx.x;
    const int iCol0T = ctColPerThread * iThreadInGrid;
    const int iCol1T = ctColPerThread * (iThreadInGrid + 1);

    EndPoint best;
    best.val = INT_MIN;
    best.p = Vec2i{ -1, -1 };

    const int blockOffset = blockRunProgress * (iBlockRun - (int)blockIdx.x) - (int)blockIdx.x;

    const int iRow0B = blockOffset;
    const int iRow1B = iRow0B + blockRunProgress;

    const int threadRowProgress = 1;
    const int threadCurRowOffset = threadRowProgress * (-(int)threadIdx.x) - (int)threadIdx.x;
    const int lastThreadCurRowOffset = threadRowProgress * (-((int)blockDim.x - 1)) - ((int)blockDim.x - 1);

    const int iRow0T = iRow0B + threadCurRowOffset;
    const int iRow1T = iRow1B + threadCurRowOffset - lastThreadCurRowOffset;

    const int iRow0S = max(0, iRow0B);
    const int iRow1S = min(ctRow, iRow1B);
    
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
        if (iRow0S <= iRow && iRow < iRow1S) {
            int prevColRow = col[iRow - 1];
            // TODO: fix this
            const bool end = iRow == iRow1S - 1 || iRow == iRow1S - 1;
            int prevCol = end && iThreadInGrid != 0 ? row[iCol0T - 1] : col[iRow];
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

                if (r > best.val && iCol < ctColReal) {
                    best.val = r;
                    best.p = Vec2i{ iRow, iCol };
                }
            }
            col[iRow - 1] = prevColRow;
            valBestInRow[iRow] = valBestInRowLocal;
            if (end) {
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

    const int ctColR = ctBlocks * ctThreadsPerBlock * ctColPerThread - 1;

    const int iLastBlock = ctBlocks - 1;
    const int ctBlockRuns = ((ctRow + iLastBlock) + blockRunProgress - 1) / blockRunProgress + iLastBlock;

    for (int iBlockRun = 0; iBlockRun < ctBlockRuns; ++iBlockRun) {
        kernel<<<ctBlocks, ctThreadsPerBlock, ctThreadPerFullBlock * ctColPerThread  * sizeof(LocalRowData)>>>(
            iBlockRun, blockRunProgress, ctColPerThread, ctRow, ctColR, ctCol, defVal,
            ca->dRow, ca->dCol, ca->dValBestInRow, ca->dValBestInCol, x, y, ca->dBest, sc);

        CUDA_CHECK(cudaGetLastError());

        CUDA_CHECK(cudaDeviceSynchronize());
    }

    return ctTotalThreads;
}
