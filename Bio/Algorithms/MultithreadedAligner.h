#pragma once

#include <Base/Core.h>

#include <algorithm>
#include <iostream>
#include <cassert>
#include <thread>
#include <atomic>
#include <functional>

class MultithreadedAligner
{
public:
  int ctMaxThreads;

  struct ThreadData {
    struct RowData {
      int row;
      int valBestInCol;
      char y;
    };

    vector<ThreadData> tds;

    vector<thread> ts;

    EndPoint _alignMultithreaded(bool rev, int iRow0, int iRow1, int iCol0, int iCol1, int rowStart, int *row,
        vector<int> &valBestInRow, vector<int> &valBestInCol, const vector<char> &x, const vector<char> &y, const Scoring &sc)
    {
      const int dir = rev ? -1 : 1;
      tds[0].col[iRow0 - dir] = row[iCol0 - dir];
      for (int iRow = iRow0; iRow != iRow1 + dir; iRow += dir) {
        tds[0].col[iRow] = rowStart;
      }
      tds[0].iRow = iRow1;

      for (int iThread = 1; iThread <= ctThreads; ++iThread) {
        tds[iThread].iRow = iRow0 - dir;
        tds[iThread].best = EndPoint{};
      }

          const char xLocal = x[iRow];
          int valBestInRowLocal = valBestInRow[iRow];
          for (int iCol = iThCol0; iCol < iThCol1; ++iCol) {
            int r = rowStart;

            auto &lrd = localRowData[iCol];

            const int prevRow = lrd.row;

            valBestInRowLocal = max(valBestInRowLocal - scLocal.k, prevCol - scLocal.b);
            lrd.valBestInCol = max(lrd.valBestInCol - scLocal.k, prevRow - scLocal.b);

            r = max(r, valBestInRowLocal);
            r = max(r, lrd.valBestInCol);

            r = max(r, prevColRow + scLocal.match(xLocal, lrd.y));

            prevColRow = prevRow;
            lrd.row = r;
            prevCol = r;

            best.Add(r, Vec2i{ iRow, iCol });
          }
          valBestInRow[iRow] = valBestInRowLocal;
          tds[iThread].col[iRow] = localRowData[iThCol1 - 1].row;
          tds[iThread].iRow = iRow;
        }
        for (int iCol = iThCol0; iCol < iThCol1; ++iCol) {
          row[iCol] = localRowData[iCol].row;
          valBestInCol[iCol] = localRowData[iCol].valBestInCol;
        }
        r.Add(tds[iThread].best);
      }
    }

  void init(int ctRow, int ctCol)
  {
    ctMaxThreads = min(ctCol, 8);

    tds.resize(ctMaxThreads + 1);
    for (auto &td : tds) {
      td._col.resize(ctRow + 3);
      td.col = td._col.data() + 1;
      td.localRowData.resize(ctCol + 2);
    }
    ts.resize(ctMaxThreads + 1);
    for (int iThread = 1; iThread <= ctMaxThreads; ++iThread) {
      ts[iThread] = thread([iThread, this]() {
        while (!tds[iThread].die) {
          if (tds[iThread].go) {
            tds[iThread].f();
            tds[iThread].go = 0;
            } else {
            this_thread::yield();
            }
            }
            });
      }
    }

  void destroy()
  {
    for (int iThread = 1; iThread <= ctMaxThreads; ++iThread) {
      tds[iThread].die = 1;
    }
    for (int iThread = 1; iThread <= ctMaxThreads; ++iThread) {
      ts[iThread].join();
    }
};
