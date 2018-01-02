#pragma once

#include <Core.h>

class NaiveCubeAA : public AlignmentAlgorithm
{
public:
  Alignment align(const string &a, const string &b, const Scoring &sc)
  {
    vector<char> x = encode(a);
    vector<char> y = encode(b);

    int xn = a.size();
    int yn = b.size();

    Mat2i dp;
    Mat2i dpr;
    dp.resize(xn + 1, yn + 1);
    dpr.resize(xn + 1, yn + 1);
    memset(dp.t.data(), 0, dp.t.size() * sizeof(char));
    for (int i = 0; i <= xn; ++i) {
      dpr[Vec2i{ i, 0 }] = INT_MAX;
    }
    for (int j = 0; j <= yn; ++j) {
      dpr[Vec2i{ 0, j }] = INT_MAX;
    }

    for (int iRow = 1; iRow <= xn; ++iRow) {
      for (int iCol = 1; iCol <= yn; ++iCol) {
        int ri = INT_MAX;
        int r = 0;

        for (int gap = 1; gap <= iRow; ++gap) {
          int mr = dp[Vec2i{ iRow - gap, iCol }] - (sc.b + sc.k * (gap - 1));
          if (mr > r) {
            r = mr;
            ri = gap;
          }
        }
        for (int gap = 1; gap <= iCol; ++gap) {
          int mr = dp[Vec2i{ iRow, iCol - gap }] - (sc.b + sc.k * (gap - 1));
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

        dp[Vec2i{ iRow, iCol }] = r;
        dpr[Vec2i{ iRow, iCol }] = ri;
      }
    }

    int s = 0;
    Vec2i p = { 0, 0 };
    for (int i = 1; i <= xn; ++i) {
      for (int j = 1; j <= yn; ++j) {
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
      if (ri == INT_MAX) {
        break;
      }
      sol.matches.push_back(p);
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
