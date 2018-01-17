#include <Base/Core.h>

#include <vector>

class CudaAligner
{
public:
    EndPoint _alignCuda(int ctRow, int ctCol, int iX, int iY, bool rev, int rowStart, int *row,
        vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc);

    void init(int ctRow, const vector<char> &x, const vector<char> &xr,
        int ctCol, const vector<char> &y, const vector<char> &yr);

    void destroy();

    // host memory
    vector<int> _col;
    int *col;
    vector<EndPoint> bests;

    // device memory
    int *_dRow;
    int *dRow;
    int *_dCol;
    int *dCol;
    char *_dX;
    char *dX;
    char *_dXr;
    char *dXr;
    char *_dY;
    char *dY;
    char *_dYr;
    char *dYr;
    int *dValBestInRow;
    int *_dValBestInCol;
    int *dValBestInCol;
    EndPoint *dBest;

};
