#include <Base/Core.h>

#include <vector>

class CudaAligner
{
public:
    EndPoint _alignCuda(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
        vector<int> &valBestInRow, vector<int> &valBestInCol, const vector<char> &x, const vector<char> &y, const Scoring &sc);

    void init(int ctRow, const vector<char> &x, int ctCol, const vector<char> &y);

    void destroy();

    // host memory
    vector<int> _col;
    int *col;
    vector<EndPoint> bests;

    // device memory
    int *dRow;
    int *dCol;
    char *dX;
    char *dY;
    int *dValBestInRow;
    int *dValBestInCol;
    EndPoint *dBest;

};
