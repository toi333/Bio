#pragma once

#include <Base/Core.h>
#include <Bio/Algorithms/MultithreadedAligner.h>
#include <CudaAligner/CudaAligner.h>

#include <algorithm>
#include <iostream>
#include <cassert>

class HirschbergAA : public AlignmentAlgorithm
{
public:
  vector<char> x, y;
  vector<char> xr, yr;
  int ctRow, ctCol;

  vector<int> rowFwdGap;
  vector<int> _rowBak;
  int *rowBak;
  vector<int> rowBakGap;

  MultithreadedAligner mta;
  CudaAligner ca;

  EndPoint _alignSingleThread(int ctRow, int ctCol, char *x, char *y, int rowStart, int *row,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    EndPoint best;

    for (int iRow = 0; iRow < ctRow; ++iRow) {
      int prevColRow = row[-1];
      row[-1] = rowStart;
      int valBestInRowLocal = valBestInRow[iRow];
      for (int iCol = 0; iCol < ctCol; ++iCol) {
        int r = rowStart;

        valBestInRowLocal = max(valBestInRowLocal - sc.k, row[iCol - 1] - sc.b);
        valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, row[iCol] - sc.b);

        r = max(r, valBestInRowLocal);
        r = max(r, valBestInCol[iCol]);

        r = max(r, prevColRow + sc.match(x[iRow], y[iCol]));

        prevColRow = row[iCol];
        row[iCol] = r;

        best.Add(r, Vec2i{ iRow, iCol });
      }
      valBestInRow[iRow] = valBestInRowLocal;
    }

    return best;
  }

  EndPoint _align(int ctRow, int ctCol, bool rev, char *x, char *y, int rowStart, int *row,
      vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
#if 0
    vector<int> rowCopy;
    for (int iCol = -1; iCol < ctCol; ++iCol) {
        rowCopy.push_back(row[iCol]);
    }
    vector<int> valBestInColCopy;
    for (int iCol = 0; iCol < ctCol; ++iCol) {
        valBestInColCopy.push_back(valBestInCol[iCol]);
    }
    vector<int> valBestInRowCopy;
    for (int iRow = 0; iRow < ctRow; ++iRow) {
        valBestInRowCopy.push_back(valBestInRow[iRow]);
    }

    EndPoint caBest = ca._alignCuda(ctRow, ctCol, x - (rev ? xr.data() : this->x.data()), y - (rev ? yr.data() : this->y.data()),
        rev, rowStart, row, valBestInRow, valBestInCol, sc);
    //EndPoint caBest = mta._alignMultithreaded(ctRow, ctCol, x, y, rowStart, row, valBestInRow, valBestInCol, sc);

    //return caBest;

    vector<int> rowCuda;
    for (int iCol = -1; iCol < ctCol; ++iCol) {
        rowCuda.push_back(row[iCol]);
    }
    vector<int> valBestInColCuda;
    for (int iCol = 0; iCol < ctCol; ++iCol) {
        valBestInColCuda.push_back(valBestInCol[iCol]);
    }

    for (int iCol = -1; iCol < ctCol; ++iCol) {
        row[iCol] = rowCopy[iCol+1];
    }
    for (int iCol = 0; iCol < ctCol; ++iCol) {
        valBestInCol[iCol] = valBestInColCopy[iCol];
    }
    for (int iRow = 0; iRow < ctRow; ++iRow) {
        valBestInRow[iRow] = valBestInRowCopy[iRow];
    }

    EndPoint stBest = _alignSingleThread(ctRow, ctCol, x, y, rowStart, row, valBestInRow, valBestInCol, sc);

    for (int iCol = -1; iCol < ctCol; ++iCol) {
        assert(row[iCol] == rowCuda[iCol+1]);
    }
    for (int iCol = 0; iCol < ctCol; ++iCol) {
        assert(valBestInCol[iCol] == valBestInColCuda[iCol]);
    }

    assert(stBest.val == caBest.val);
    assert(stBest.p.x == caBest.p.x);
    assert(stBest.p.y == caBest.p.y);

    return stBest;
#else
#if 0
      return ca._alignCuda(ctRow, ctCol, x - (rev ? xr.data() : this->x.data()), y - (rev ? yr.data() : this->y.data()),
          rev, rowStart, row, valBestInRow, valBestInCol, sc);
#else
    if (ctRow < 32 || ctCol < 32) {
      return _alignSingleThread(ctRow, ctCol, x, y, rowStart, row, valBestInRow, valBestInCol, sc);
    }
    #if 1
        if (ctRow < 2000 || ctCol < 2000) {
            return mta._alignMultithreaded(ctRow, ctCol, x, y, rowStart, row, valBestInRow, valBestInCol, sc);
        }
        return ca._alignCuda(ctRow, ctCol, x - (rev ? xr.data() : this->x.data()), y - (rev ? yr.data() : this->y.data()),
            rev, rowStart, row, valBestInRow, valBestInCol, sc);
    #else
        return mta._alignMultithreaded(ctRow, ctCol, x, y, rowStart, row, valBestInRow, valBestInCol, sc);
    #endif
#endif
#endif
  }

  EndPoint alignLastRow(int ctCol, char *x, char *y, int rowStart, int *row, vector<int> &rowGap,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    EndPoint best;

    int prevColRow = rowStart;
    for (int iCol = 0; iCol < ctCol; ++iCol) {
      valBestInRow[0] = max(valBestInRow[0] - sc.k, row[iCol - 1] - sc.b);
      valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, row[iCol] - sc.b);

      int rGap = valBestInCol[iCol];
      int rNoGap = max(prevColRow + sc.match(x[0], y[iCol]), valBestInRow[0]);
      int r = max(max(rGap, rNoGap), rowStart);

      prevColRow = row[iCol];

      rowGap[iCol] = rGap;
      row[iCol] = r;

      best.Add(r, Vec2i{ 0, iCol });
    }

    return best;
  }

  void alignRec(int iRow0, int iRow1, int iCol0, int iCol1, bool gapStart, bool gapEnd,
    int *row, vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc, Alignment &alig)
  {
    assert(iRow0 <= iRow1);
    assert(iCol0 <= iCol1);

    const int negInf = -(1 << 28);

    const bool startFree = iRow0 == 0 && iCol0 == 0;
    const bool endFree = iRow1 == ctRow && iCol1 == ctCol;

    if (iRow0 == iRow1 || iCol0 == iCol1) {
      if (startFree) {
        iRow0 = iRow1;
        iCol0 = iCol1;
      }
      if (endFree) {
        iRow1 = iRow0;
        iCol1 = iCol0;
      }
      if (iCol0 == iCol1) {
        for (int iRow = iRow0; iRow <= iRow1; ++iRow) {
          alig.matches.push_back(Vec2i{ iRow, iCol0 });
        }
      } else {
        for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
          alig.matches.push_back(Vec2i{ iRow0, iCol });
        }
      }
      return;
    }

    //if (iRow0 + 1 == iRow1) {
    //  const int gapScore = -((iCol1 - iCol0) * sc.k + sc.b + sc.k + (gapStart ? 0 : sc.b));
    //  int bestMatchScore = negInf;
    //  int bestCol = -1;
    //  for (int iCol = iCol0 + 1; iCol <= iCol1; ++iCol) {
    //    const int matchScore = sc.m[x[iRow1]][y[iCol]];
    //    if (matchScore > bestMatchScore) {
    //      bestMatchScore = matchScore;
    //      bestCol = iCol;
    //    }
    //  }
    //  bestMatchScore -= (iCol1 - iCol0 - 1) * sc.k + 2 * sc.b;
    //  if (bestMatchScore > gapScore) {
    //    alignRec(iRow0, iRow1, iCol0, iCol0, gapStart, row, valBestInRow, valBestInCol, sc, alig);
    //    alignRec(iRow1, iRow1, iCol0, iCol1, gapStart, row, valBestInRow, valBestInCol, sc, alig);
    //  } else {
    //    alignRec(iRow0, iRow1, iCol0, iCol0, gapStart, row, valBestInRow, valBestInCol, sc, alig);
    //    alignRec(iRow1, iRow1, iCol0, iCol1, gapStart, row, valBestInRow, valBestInCol, sc, alig);
    //  }
    //  return;
    //}

    int iMidRow = (iRow0 + iRow1) / 2;

    // First half forward alignment
    const int startVal = startFree ? 0 : negInf;
    row[-1] = startVal;
    row[0] = 0;
    for (int iCol = iCol0 + 1; iCol <= iCol1; ++iCol) {
      row[iCol-iCol0] = max(-(sc.b + (iCol - iCol0 - 1) * sc.k), startVal);
    }
    for (int iRow = iRow0; iRow <= iMidRow; ++iRow) {
      valBestInRow[iRow - iRow0] = negInf;
    }
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      valBestInCol[iCol - iCol0] = negInf;
    }
    if (gapStart) {
      valBestInCol[0] = 0;
    }
    EndPoint bestFwd{ 0, { iRow0, iCol0} };
    if (iRow0 + 1 <= iMidRow - 1) {
      EndPoint b = _align(iMidRow - iRow0 - 1, iCol1 - iCol0 + 1, false,
          x.data() + iRow0 + 1, y.data() + iCol0, startVal, row, valBestInRow, valBestInCol, sc);
      b.p.x += iRow0 + 1;
      b.p.y += iCol0;
      bestFwd.Add(b);
    }

    if (iRow0 != iMidRow) {
      valBestInRow[0] = negInf;
      EndPoint b = alignLastRow(iCol1 - iCol0 + 1, x.data() + iMidRow, y.data() + iCol0,
          startVal, row, rowFwdGap, valBestInRow, valBestInCol, sc);
      b.p.x += iMidRow;
      b.p.y += iCol0;
      bestFwd.Add(b);
    } else {
      for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
        rowFwdGap[iCol - iCol0] = negInf;
      }
    }


    // Last half backward alignment
    const int endVal = endFree ? 0 : negInf;
    rowBak[-1] = endVal;
    rowBak[0] = 0;
    for (int iCol = iCol0; iCol < iCol1; ++iCol) {
      rowBak[iCol1-iCol] = max(-(sc.b + (iCol1 - iCol - 1) * sc.k), endVal);
    }
    for (int iRow = iMidRow; iRow <= iRow1; ++iRow) {
      valBestInRow[iRow-iMidRow] = negInf;
    }
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      valBestInCol[iCol-iCol0] = negInf;
    }
    if (gapEnd) {
      valBestInCol[0] = 0;
    }
    EndPoint bestBak{ 0, { iRow1, iCol1 } };
    if (iMidRow + 1 <= iRow1 - 1) {
      EndPoint b = _align(iRow1 - iMidRow - 1, iCol1 - iCol0 + 1, true,
          xr.data() + ctRow - (iRow1 - 1), yr.data() + ctCol - iCol1, endVal, rowBak, valBestInRow, valBestInCol, sc);
      b.p.x = iRow1 - 1 - b.p.x;
      b.p.y = iCol1 - b.p.y;
      bestBak.Add(b);
    }

    if (iRow1 != iMidRow) {
      valBestInRow[0] = negInf;
      EndPoint b = alignLastRow(iCol1 - iCol0 + 1, xr.data() + ctRow - iMidRow, yr.data() + ctCol - iCol1,
          endVal, rowBak, rowBakGap, valBestInRow, valBestInCol, sc);
      b.p.x = iMidRow - b.p.x;
      b.p.y = iCol1 - b.p.y;
      bestBak.Add(b);
    } else {
      for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
        rowBakGap[iCol - iCol0] = negInf;
      }
    }


    // Merge two halves
    int iColMax = -1;
    int valMax = negInf;
    bool gap = false;
    bool gapDown = false;
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      const int valGap = rowFwdGap[iCol - iCol0] + rowBakGap[iCol1 - iCol] + sc.b - sc.k;
      if (valGap >= valMax) {
        valMax = valGap;
        iColMax = iCol;
        gap = true;
      }
      const int valNoGap = row[iCol - iCol0] + rowBak[iCol1 - iCol];
      if (valNoGap >= valMax) {
        valMax = valNoGap;
        iColMax = iCol;
        gap = false;
        gapDown = rowBakGap[iCol1 - iCol] == rowBak[iCol1 - iCol];
      }
    }


    // If last half should be discarded
    if (endFree && bestFwd.val >= valMax && (!startFree || bestFwd.val >= bestBak.val)) {
      if (startFree) {
        alig.score = bestFwd.val;
      }
      if (bestFwd.val > 0) {
        alignRec(iRow0, bestFwd.p.x - 1, iCol0, bestFwd.p.y - 1, gapStart, false, row, valBestInRow, valBestInCol, sc, alig);
      } else {
        assert(bestFwd.p.x == iRow0 && bestFwd.p.y == iCol0);
      }
      alig.matches.push_back(bestFwd.p);
      return;
    }

    // If first half should be discarded
    if (startFree && bestBak.val >= valMax && (!endFree || bestBak.val >= bestFwd.val)) {
      if (endFree) {
        alig.score = bestBak.val;
      }
      alig.matches.push_back(bestBak.p);
      if (bestBak.val > 0) {
        alignRec(bestBak.p.x + 1, iRow1, bestBak.p.y + 1, iCol1, false, gapEnd, row, valBestInRow, valBestInCol, sc, alig);
      } else {
        assert(bestBak.p.x == iRow1 && bestBak.p.y == iCol1);
      }
      return;
    }

    // If this is the root call
    if (startFree && endFree) {
      alig.score = valMax;
    }

    // Reconstruct recursively
    if (gap) {
      alignRec(iRow0, iMidRow - 1, iCol0, iColMax, gapStart, true, row, valBestInRow, valBestInCol, sc, alig);
      alig.matches.emplace_back(Vec2i{ iMidRow, iColMax });
      alignRec(iMidRow + 1, iRow1, iColMax, iCol1, true, gapEnd, row, valBestInRow, valBestInCol, sc, alig);
    } else {
      alignRec(iRow0, iMidRow, iCol0, iColMax, gapStart, gapDown, row, valBestInRow, valBestInCol, sc, alig);
      alignRec(iMidRow + 1, iRow1, iColMax + !gapDown, iCol1, gapDown, gapEnd, row, valBestInRow, valBestInCol, sc, alig);
    }
  }

  Alignment align(const string &a, const string &b, const Scoring &sc)
  {
    x = encode(a);
    y = encode(b);

    xr = x;
    reverse(xr.begin() + 1, xr.end());
    yr = y;
    reverse(yr.begin() + 1, yr.end());

    ctRow = (int)a.size();
    ctCol = (int)b.size();

    vector<int> row(ctCol + 3, 0);

    const int negInf = -(1 << 28);

    vector<int> valBestInRow(ctRow + 1, negInf);
    vector<int> valBestInCol(ctCol + 1, negInf);

    rowFwdGap.resize(ctCol + 1);
    _rowBak.resize(ctCol + 2);
    rowBak = _rowBak.data() + 1;
    rowBakGap.resize(ctCol + 2);

    mta.init(ctRow, ctCol);
    ca.init(ctRow, x, xr, ctCol, y, yr);

    Alignment sol;
    alignRec(0, ctRow, 0, ctCol, false, false, row.data() + 1, valBestInRow, valBestInCol, sc, sol);
    sol.compress();

    mta.destroy();
    ca.destroy();

    return sol;
  }
};
