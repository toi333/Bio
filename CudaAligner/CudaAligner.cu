#include "cuda_runtime.h"
#include "device_launch_parameters.h"

#include <Base/Core.h>

#include <vector>

class CudaAligner
{
public:

    EndPoint _alignCuda(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
        vector<int> &valBestInRow, vector<int> &valBestInCol, const vector<char> &x, const vector<char> &y, const Scoring &sc)
    {
    }


    void init(int ctRow, int ctCol)
    {
    }


    void destroy()
    {
    }
};
