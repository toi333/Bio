#include <CudaAligner/CudaAligner.h>

#include <CudaAligner/Util.h>

#include <cuda_runtime.h>

#include <iostream>
#include <ctime>

extern int callKernel(int ctRow, int ctCol, char *x, char *y, int defVal, CudaAligner *ca, const Scoring &sc);

EndPoint CudaAligner::_alignCuda(int ctRow, int ctCol, int iX, int iY, bool rev, int rowStart, int *row,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
//EndPoint CudaAligner::_alignCuda(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
//    vector<int> &valBestInRow, vector<int> &valBestInCol, const vector<char> &x, const vector<char> &y, const Scoring &sc)
{
    col[-1] = row[-1];
    for (int iRow = 0; iRow < ctRow; ++iRow) {
        col[iRow] = rowStart;
    }

    CUDA_CHECK(cudaMemcpy(dRow, row, ctCol * sizeof(*row), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dCol - 1, col - 1, (ctRow + 1) * sizeof(*col), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dValBestInCol, valBestInCol.data(), ctCol * sizeof(*dValBestInCol), cudaMemcpyHostToDevice));
    CUDA_CHECK(cudaMemcpy(dValBestInRow, valBestInRow.data(), ctRow * sizeof(*dValBestInRow), cudaMemcpyHostToDevice));

    //cout << "kernel " << iCol0 << ' ' << iCol1 << endl;
    static int totalTime = 0;
    const int start = clock();

    int ctBest = callKernel(ctRow, ctCol, (rev ? _dXr : _dX) + iX, (rev ? _dYr : _dY) + iY, rowStart, this, sc);

    CUDA_CHECK(cudaGetLastError());

    CUDA_CHECK(cudaDeviceSynchronize());

    const int tm = clock() - start;
    totalTime += tm;
    if (tm > 1000) {
        cout << "time: " << (float)tm / CLOCKS_PER_SEC << ' ' << (float)totalTime / CLOCKS_PER_SEC << endl;
    }

    CUDA_CHECK(cudaMemcpy(row, dRow, ctCol * sizeof(*row), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(valBestInCol.data(), dValBestInCol, ctCol * sizeof(*dValBestInCol), cudaMemcpyDeviceToHost));

    CUDA_CHECK(cudaMemcpy(bests.data(), dBest, ctBest * sizeof(*dBest), cudaMemcpyDeviceToHost));
    EndPoint best = bests[0];
    for (int i = 1; i < ctBest; ++i) {
        const EndPoint &b = bests[i];
        if (b.val > best.val || b.val == best.val &&
            (b.p.x < best.p.x || b.p.x == best.p.x && b.p.y < best.p.y)) {
            best = b;
        }
    }

    return best;
}


void CudaAligner::init(int ctRow, const vector<char> &x, const vector<char> &xr,
    int ctCol, const vector<char> &y, const vector<char> &yr)
{
    // TODO: WTF
    const int ctMaxBests = 1 << 18;
    const int ctMaxBlocks = (ctCol + (1 << 7) - 1) / (1 << 7) + 100;

    _col.resize(ctRow + 3);
    col = _col.data() + 1;

    bests.resize(ctMaxBests);

    // TODO: check sentinels?
    // TODO: allocations alignment?
    cudaMalloc(&_dRow, (ctCol + 2 + 2 * ctMaxBlocks) * sizeof(*dRow));
    dRow = _dRow + ctMaxBlocks;
    cudaMalloc(&_dCol, (ctRow + 2) * sizeof(*dCol));
    dCol = _dCol + 1;
    cudaMalloc(&_dX, (ctRow + 1) * sizeof(*dX));
    cudaMalloc(&_dY, (ctCol + 1 + 2 * ctMaxBlocks) * sizeof(*dY));
    cudaMalloc(&_dXr, (ctRow + 1) * sizeof(*dX));
    cudaMalloc(&_dYr, (ctCol + 1 + 2 * ctMaxBlocks) * sizeof(*dY));
    dX = _dX + 1;
    dXr = _dXr + 1;
    dY = _dY + 1 + ctMaxBlocks;
    dYr = _dYr + 1 + ctMaxBlocks;
    cudaMalloc(&dValBestInRow, (ctRow + 2) * sizeof(*dValBestInRow));
    cudaMalloc(&_dValBestInCol, (ctCol + 2 + 2 * ctMaxBlocks) * sizeof(*dValBestInCol));
    dValBestInCol = _dValBestInCol + ctMaxBlocks;
    cudaMalloc(&dBest, ctMaxBests * sizeof(*dBest));

    cudaMemcpy(_dX, x.data(), (ctRow + 1) * sizeof(*dX), cudaMemcpyHostToDevice);
    cudaMemcpy(_dY, y.data(), (ctCol + 1) * sizeof(*dY), cudaMemcpyHostToDevice);
    cudaMemcpy(_dXr, xr.data(), (ctRow + 1) * sizeof(*dX), cudaMemcpyHostToDevice);
    cudaMemcpy(_dYr, yr.data(), (ctCol + 1) * sizeof(*dY), cudaMemcpyHostToDevice);
}

void CudaAligner::destroy()
{
    cudaFree(dBest);
    cudaFree(_dValBestInCol);
    cudaFree(dValBestInRow);
    cudaFree(_dY);
    cudaFree(_dX);
    cudaFree(_dYr);
    cudaFree(_dXr);
    cudaFree(_dCol);
    cudaFree(_dRow);
}
