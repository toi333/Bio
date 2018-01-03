#pragma once

#include <Core.h>

class SquareAA : public AlignmentAlgorithm
{
public:
  Alignment align(const string &a, const string &b, const Scoring &sc)
  {
    vector<char> x = encode(a);
    vector<char> y = encode(b);

    int ctRow = (int)a.size();
    int ctCol = (int)b.size();

    Mat2i dp;
    Mat2i dpr;
    dp.resize(ctRow + 1, ctCol + 1);
    dpr.resize(ctRow + 1, ctCol + 1);
    memset(dp.t.data(), 0, dp.t.size() * sizeof(char));
    for (int i = 0; i <= ctRow; ++i) {
      dpr[Vec2i{ i, 0 }] = INT_MAX;
    }
    for (int j = 0; j <= ctCol; ++j) {
      dpr[Vec2i{ 0, j }] = INT_MAX;
    }

    vector<int> idxBestInRow(ctRow + 1, 0);
    vector<int> idxBestInCol(ctCol + 1, 0);

    for (int iRow = 1; iRow <= ctRow; ++iRow) {
      for (int iCol = 1; iCol <= ctCol; ++iCol) {
        int ri = INT_MAX;
        int r = 0;

        {
          int gap = iRow - idxBestInCol[iCol];
          int mr = dp[Vec2i{ idxBestInCol[iCol], iCol }] - (sc.b + sc.k * (gap - 1));
          if (mr > r) {
            r = mr;
            ri = gap;
          }
        }
        {
          int gap = iCol - idxBestInRow[iRow];
          int mr = dp[Vec2i{ iRow, idxBestInRow[iRow] }] - (sc.b + sc.k * (gap - 1));
          if (mr > r) {
            r = mr;
            ri = -gap;
          }
        }

        int mr = dp[Vec2i{ iRow - 1, iCol - 1 }] + sc.m[x[iRow]][y[iCol]];
        if (mr > r) {
          r = mr;
          ri = 0;
        }

        if (dp[Vec2i{ iRow, idxBestInRow[iRow] }] + sc.k * idxBestInRow[iRow] < r + sc.k * iCol) {
          idxBestInRow[iRow] = iCol;
        }
        if (dp[Vec2i{ idxBestInCol[iCol], iCol }] + sc.k * idxBestInCol[iCol] < r + sc.k * iRow) {
          idxBestInCol[iCol] = iRow;
        }
        dp[Vec2i{ iRow, iCol }] = r;
        dpr[Vec2i{ iRow, iCol }] = ri;
      }
    }

    int s = 0;
    Vec2i p = { 0, 0 };
    for (int i = 1; i <= ctRow; ++i) {
      for (int j = 1; j <= ctCol; ++j) {
        if (dp[Vec2i{ i, j }] > s) {
          s = dp[Vec2i{ i, j }];
          p = Vec2i{ i, j };
        }
      }
    }

    if (s == 0) {
      return Alignment();
    }

    Alignment sol;
    sol.score = s;
    while (1337) {
      int ri = dpr[p];
      sol.matches.push_back(p);
      if (ri == INT_MAX) {
        break;
      }
      if (ri == 0) {
        p = { p.x - 1, p.y - 1 };
      } else if (ri > 0) {
        p = { p.x - ri, p.y };
      } else {
        p = { p.x, p.y + ri };
      }
    }
    reverse(sol.matches.begin(), sol.matches.end());
    return sol;
  }
};
