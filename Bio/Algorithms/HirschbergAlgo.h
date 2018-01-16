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
  int ctRow, ctCol;

  vector<int> rowFwdGap;
  vector<int> rowBak;
  vector<int> rowBakGap;

  MultithreadedAligner mta;
  CudaAligner ca;

  EndPoint _alignSingleThread(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    const int dir = rev ? -1 : 1;

    iRow1 += dir;
    iCol1 += dir;

    EndPoint best;

    for (int iRow = iRow0; iRow != iRow1; iRow += dir) {
      int prevColRow = row[iCol0 - dir];
      row[iCol0 - dir] = rowStart;
      int valBestInRowLocal = valBestInRow[iRow];
      for (int iCol = iCol0; iCol != iCol1; iCol += dir) {
        int r = rowStart;

        valBestInRowLocal = max(valBestInRowLocal - sc.k, row[iCol - dir] - sc.b);
        valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, row[iCol] - sc.b);

        r = max(r, valBestInRowLocal);
        r = max(r, valBestInCol[iCol]);

        r = max(r, prevColRow + sc.match(x[iRow + rev], y[iCol + rev]));

        prevColRow = row[iCol];
        row[iCol] = r;

        best.Add(r, Vec2i{ iRow, iCol });
      }
      valBestInRow[iRow] = valBestInRowLocal;
    }

    return best;
  }

  EndPoint _align(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    //int dir = rev ? -1 : 1;

    //vector<int> rowCopy;
    //for (int iCol = iCol0 - dir; iCol != iCol1 + dir; iCol += dir) {
    //    rowCopy.push_back(row[iCol]);
    //}
    //vector<int> valBestInColCopy;
    //for (int iCol = iCol0; iCol != iCol1 + dir; iCol += dir) {
    //    valBestInColCopy.push_back(valBestInCol[iCol]);
    //}

    //EndPoint caBest = ca._alignCuda(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, x, y, sc);

    ////return caBest;

    //vector<int> rowCuda;
    //for (int iCol = iCol0 - dir; iCol != iCol1 + dir; iCol += dir) {
    //    rowCuda.push_back(row[iCol]);
    //}
    //vector<int> valBestInColCuda;
    //for (int iCol = iCol0; iCol != iCol1 + dir; iCol += dir) {
    //    valBestInColCuda.push_back(valBestInCol[iCol]);
    //}

    //for (int iCol = iCol0 - dir, j = 0; iCol != iCol1 + dir; iCol += dir, ++j) {
    //    row[iCol] = rowCopy[j];
    //}
    //for (int iCol = iCol0, j = 0; iCol != iCol1 + dir; iCol += dir, ++j) {
    //    valBestInCol[iCol] = valBestInColCopy[j];
    //}

    //EndPoint stBest = _alignSingleThread(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, sc);

    //for (int iCol = iCol0 - dir, j = 0; iCol != iCol1 + dir; iCol += dir, ++j) {
    //    assert(row[iCol] == rowCuda[j]);
    //}
    //for (int iCol = iCol0, j = 0; iCol != iCol1 + dir; iCol += dir, ++j) {
    //    assert(valBestInCol[iCol] == valBestInColCuda[j]);
    //}

    //assert(stBest.val == caBest.val);
    //assert(stBest.p.x == caBest.p.x);
    //assert(stBest.p.y == caBest.p.y);

    //return stBest;
#if 0
    return _alignSingleThread(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, sc);
#else
    if (abs(iRow0 - iRow1) < 32 || abs(iCol0 - iCol1) < 32) {
      return _alignSingleThread(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, sc);
    }
    #if 1
        if (abs(iRow0 - iRow1) < 2000 || abs(iCol0 - iCol1) < 2000) {
            return mta._alignMultithreaded(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, x, y, sc);
        }
        return ca._alignCuda(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, x, y, sc);
    #else
        return mta._alignMultithreaded(rev, iRow0, iRow1, iCol0, iCol1, rowStart, row, valBestInRow, valBestInCol, x, y, sc);
    #endif
#endif
  }

  EndPoint alignLastRow(int iRow, int iCol0, int iCol1, int rowStart, int *row, vector<int> &rowGap,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    const int colRev = iCol0 < iCol1 ? 0 : 1;
    const int rowRev = colRev;

    const int colDir = !colRev ? 1 : -1;

    iCol1 += colDir;

    EndPoint best;

    int prevColRow = rowStart;
    for (int iCol = iCol0; iCol != iCol1; iCol += colDir) {
      valBestInRow[iRow] = max(valBestInRow[iRow] - sc.k, row[iCol - colDir] - sc.b);
      valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, row[iCol] - sc.b);

      int rGap = valBestInCol[iCol];
      int rNoGap = max(prevColRow + sc.match(x[iRow + rowRev], y[iCol + colRev]), valBestInRow[iRow]);
      int r = max(max(rGap, rNoGap), rowStart);

      prevColRow = row[iCol];

      rowGap[iCol] = rGap;
      row[iCol] = r;

      best.Add(r, Vec2i{ iRow, iCol });
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
    row[iCol0 - 1] = startVal;
    row[iCol0] = 0;
    for (int iCol = iCol0 + 1; iCol <= iCol1; ++iCol) {
      row[iCol] = max(-(sc.b + (iCol - iCol0 - 1) * sc.k), startVal);
    }
    for (int iRow = iRow0; iRow <= iMidRow; ++iRow) {
      valBestInRow[iRow] = negInf;
    }
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      valBestInCol[iCol] = negInf;
    }
    if (gapStart) {
      valBestInCol[iCol0] = 0;
    }
    EndPoint bestFwd{ 0, { iRow0, iCol0} };
    if (iRow0 + 1 <= iMidRow - 1) {
      bestFwd.Add(_align(false, iRow0 + 1, iMidRow - 1, iCol0, iCol1, startVal, row, valBestInRow, valBestInCol, sc));
    }

    if (iRow0 != iMidRow) {
      bestFwd.Add(alignLastRow(iMidRow, iCol0, iCol1, startVal, row, rowFwdGap, valBestInRow, valBestInCol, sc));
    } else {
      for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
        rowFwdGap[iCol] = negInf;
      }
    }


    // Last half backward alignment
    const int endVal = endFree ? 0 : negInf;
    for (int iCol = iCol0; iCol < iCol1; ++iCol) {
      rowBak[iCol] = max(-(sc.b + (iCol1 - iCol - 1) * sc.k), endVal);
    }
    rowBak[iCol1] = 0;
    rowBak[iCol1 + 1] = endVal;
    for (int iRow = iMidRow; iRow <= iRow1; ++iRow) {
      valBestInRow[iRow] = negInf;
    }
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      valBestInCol[iCol] = negInf;
    }
    if (gapEnd) {
      valBestInCol[iCol1] = 0;
    }
    EndPoint bestBak{ 0, { iRow1, iCol1 } };
    if (iMidRow + 1 <= iRow1 - 1) {
      bestBak.Add(_align(true, iRow1 - 1, iMidRow + 1, iCol1, iCol0, endVal, rowBak.data(), valBestInRow, valBestInCol, sc));
    }

    if (iRow1 != iMidRow) {
      bestBak.Add(alignLastRow(iMidRow, iCol1, iCol0, endVal, rowBak.data(), rowBakGap, valBestInRow, valBestInCol, sc));
    } else {
      for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
        rowBakGap[iCol] = negInf;
      }
    }


    // Merge two halves
    int iColMax = -1;
    int valMax = negInf;
    bool gap = false;
    bool gapDown = false;
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      const int valGap = rowFwdGap[iCol] + rowBakGap[iCol] + sc.b - sc.k;
      if (valGap >= valMax) {
        valMax = valGap;
        iColMax = iCol;
        gap = true;
      }
      const int valNoGap = row[iCol] + rowBak[iCol];
      if (valNoGap >= valMax) {
        valMax = valNoGap;
        iColMax = iCol;
        gap = false;
        gapDown = rowBakGap[iCol] == rowBak[iCol];
      }
    }


    // If last half should be discarded
    if (endFree && bestFwd.val >= valMax && (!startFree || bestFwd.val >= bestBak.val)) {
      if (startFree) {
        alig.score = bestFwd.val;
      }
      if (bestFwd.val > 0) {
        alignRec(iRow0, bestFwd.p.x - 1, iCol0, bestFwd.p.y - 1, false, false, row, valBestInRow, valBestInCol, sc, alig);
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
        alignRec(bestBak.p.x + 1, iRow1, bestBak.p.y + 1, iCol1, false, false, row, valBestInRow, valBestInCol, sc, alig);
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

    ctRow = (int)a.size();
    ctCol = (int)b.size();

    vector<int> row(ctCol + 3, 0);

    const int negInf = -(1 << 28);

    vector<int> valBestInRow(ctRow + 1, negInf);
    vector<int> valBestInCol(ctCol + 1, negInf);

    rowFwdGap.resize(ctCol + 1);
    rowBak.resize(ctCol + 2);
    rowBakGap.resize(ctCol + 2);

    mta.init(ctRow, ctCol);
    ca.init(ctRow, x, ctCol, y);

    Alignment sol;
    alignRec(0, ctRow, 0, ctCol, false, false, row.data() + 1, valBestInRow, valBestInCol, sc, sol);
    sol.compress();

    mta.destroy();
    ca.destroy();

    return sol;
  }
};
