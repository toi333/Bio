#pragma once

#include <Core.h>

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

  struct EndPoint
  {
    int val = -1;
    Vec2i p;

    bool Add(int val, Vec2i p)
    {
      return Add(EndPoint{ val, p });
    }

    bool Add(EndPoint ep)
    {
      if (ep.val > val) {
        val = ep.val;
        p = ep.p;
        return true;
      }
      return false;
    }
  };

  // TODO: local instead of rowStart?
  EndPoint _align(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row, vector<int> &valBestInRow, vector<int> &valBestInCol,
    const Scoring &sc)
  {
    const int dir = rev ? -1 : 1;

    iRow1 += dir;
    iCol1 += dir;

    EndPoint best;

    for (int iRow = iRow0; iRow != iRow1; iRow += dir) {
      int prevColRow = row[iCol0 - dir];
      row[iCol0 - dir] = rowStart;
      for (int iCol = iCol0; iCol != iCol1; iCol += dir) {
        int r = rowStart;

        valBestInRow[iRow] = max(valBestInRow[iRow] - sc.k, row[iCol - dir] - sc.b);
        valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, row[iCol] - sc.b);

        r = max(r, valBestInRow[iRow]);
        r = max(r, valBestInCol[iCol]);

        r = max(r, prevColRow + sc.m[x[iRow + rev]][y[iCol + rev]]);

        prevColRow = row[iCol];
        row[iCol] = r;

        best.Add(r, Vec2i{ iRow, iCol });
      }
    }

    return best;
  }
  
  pair<EndPoint, EndPoint> findEndPoints(int *row, vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    EndPoint ep1 = _align(false, 1, ctRow, 1, ctCol, 0, row, valBestInRow, valBestInCol, sc);
    if (ep1.val <= 0) {
      return make_pair(EndPoint{ 0, Vec2i{ 0, 0 } }, EndPoint{ 0, Vec2i{ 0, 0 } });
    }
    // TODO check
    memset(row, 0, (ctCol + 2) * sizeof(int));
    for (int &v : valBestInRow) {
      v = -(1 << 28);
    }
    for (int &v : valBestInCol) {
      v = -(1 << 28);
    }

    EndPoint ep0 = _align(true, ep1.p.x - 1, 0, ep1.p.y - 1, 0, 0, row, valBestInRow, valBestInCol, sc);
    return make_pair(ep0, ep1);
  }

  EndPoint alignLastRow(int iRow, int iCol0, int iCol1, int *row, vector<int> &rowGap,
    vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    const int colRev = iCol0 < iCol1 ? 0 : 1;
    const int rowRev = colRev;

    const int colDir = !colRev ? 1 : -1;

    iCol1 += colDir;

    EndPoint best;

    int prevColRow = -(1 << 28); //0;
    for (int iCol = iCol0; iCol != iCol1; iCol += colDir) {
      valBestInRow[iRow] = max(valBestInRow[iRow] - sc.k, row[iCol - colDir] - sc.b);
      valBestInCol[iCol] = max(valBestInCol[iCol] - sc.k, row[iCol] - sc.b);

      int rGap = valBestInCol[iCol];
      int rNoGap = max(prevColRow + sc.m[x[iRow + rowRev]][y[iCol + colRev]], valBestInRow[iRow]);
      int r = max(rGap, rNoGap);

      prevColRow = row[iCol];

      rowGap[iCol] = rGap;
      row[iCol] = r;

      best.Add(r, Vec2i{ iRow, iCol });
    }

    return best;
  }

  void alignRec(int iRow0, int iRow1, int iCol0, int iCol1, bool gapStart, bool gapEnd, int *row, vector<int> &valBestInRow,
    vector<int> &valBestInCol, const Scoring &sc, Alignment &alig)
  {
    assert(iRow0 <= iRow1);
    assert(iCol0 <= iCol1);

    const int negInf = -(1 << 28);

    if (iRow0 == iRow1 || iCol0 == iCol1) {
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

    row[iCol0 - 1] = negInf;
    row[iCol0] = 0;
    for (int iCol = iCol0 + 1; iCol <= iCol1; ++iCol) {
      row[iCol] = -(sc.b + (iCol - iCol0 - 1) * sc.k);
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
    if (iRow0 + 1 <= iMidRow - 1) {
      _align(false, iRow0 + 1, iMidRow - 1, iCol0, iCol1, negInf, row, valBestInRow, valBestInCol, sc);
    }

    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      rowFwdGap[iCol] = negInf;
    }

    if (iRow0 != iMidRow) {
      alignLastRow(iMidRow, iCol0, iCol1, row, rowFwdGap, valBestInRow, valBestInCol, sc);
    }


    for (int iCol = iCol0; iCol < iCol1; ++iCol) {
      rowBak[iCol] = -(sc.b + (iCol1 - iCol - 1) * sc.k);
    }
    rowBak[iCol1] = 0;
    rowBak[iCol1 + 1] = negInf;
    for (int iRow = iMidRow; iRow <= iRow1; ++iRow) {
      valBestInRow[iRow] = negInf;
    }
    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      valBestInCol[iCol] = negInf;
    }
    if (gapEnd) {
      valBestInCol[iCol1] = 0;
    }

    if (iMidRow + 1 <= iRow1 - 1) {
      _align(true, iRow1 - 1, iMidRow + 1, iCol1, iCol0, negInf, rowBak.data(), valBestInRow, valBestInCol, sc);
    }

    for (int iCol = iCol0; iCol <= iCol1; ++iCol) {
      rowBakGap[iCol] = negInf;
    }

    if (iRow1 != iMidRow) {
      alignLastRow(iMidRow, iCol1, iCol0, rowBak.data(), rowBakGap, valBestInRow, valBestInCol, sc);
    }

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

    ctRow = a.size();
    ctCol = b.size();

    vector<int> row(ctCol + 3, 0);

    const int negInf = -(1 << 28);

    vector<int> valBestInRow(ctRow + 1, negInf);
    vector<int> valBestInCol(ctCol + 1, negInf);

    auto eps = findEndPoints(row.data(), valBestInRow, valBestInCol, sc);
    const EndPoint ep0 = eps.first;
    const EndPoint ep1 = eps.second;
    assert(ep0.val == ep1.val);

    rowFwdGap.resize(ctCol + 1);
    rowBak.resize(ctCol + 2);
    rowBakGap.resize(ctCol + 2);

    Alignment sol;
    sol.score = ep0.val;
    alignRec(ep0.p.x, ep1.p.x, ep0.p.y, ep1.p.y, false, false, row.data() + 1, valBestInRow, valBestInCol, sc, sol);
    sol.compress();
    return sol;

    //EndPoint best;

    //int iMidRow = (ep0.p.x + ep1.p.x) / 2;


    //for (int &v : row) {
    //  v = negInf;
    //}
    //row[ep0.p.y] = 0;
    //for (int &v : valBestInRow) {
    //  v = negInf;
    //}
    //for (int &v : valBestInCol) {
    //  v = negInf;
    //}
    //if (iMidRow != ep0.p.x) {
    //  best.Add(_align(ep0.p.x + 1, iMidRow - 1, ep0.p.y + 1, ep1.p.y, negInf, row, valBestInRow, valBestInCol, sc));
    //}

    //vector<int> rowFwdGap(ctCol + 1, negInf);

    //best.Add(alignLastRow(iMidRow, ep0.p.y + 1, ep1.p.y, row, rowFwdGap, valBestInRow, valBestInCol, sc));


    //vector<int> rowBak(ctCol + 2, negInf);
    //rowBak[ep1.p.y] = 0;
    //for (int &v : valBestInRow) {
    //  v = negInf;
    //}
    //for (int &v : valBestInCol) {
    //  v = negInf;
    //}

    //if (ep1.p.x != iMidRow) {
    //  best.Add(_align(ep1.p.x - 1, iMidRow + 1, ep1.p.y - 1, ep0.p.y, negInf, rowBak, valBestInRow, valBestInCol, sc));
    //}

    //vector<int> rowBakGap(ctCol + 2, negInf);

    //best.Add(alignLastRow(iMidRow, ep1.p.y - 1, ep0.p.y, rowBak, rowBakGap, valBestInRow, valBestInCol, sc));


    //for (int iCol = ep0.p.y; iCol <= ep1.p.y; ++iCol) {
    //  int r = row[iCol] + rowBak[iCol];
    //  r = max(r, rowFwdGap[iCol] + rowBakGap[iCol] + sc.b - sc.k);
    //  best.Add(r, Vec2i{ iMidRow, iCol });
    //}

    //Alignment sol;
    //sol.score = best.val;
    //return sol;
  }
};
