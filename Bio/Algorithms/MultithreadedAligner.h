#pragma once

#include <Core.h>

#include <algorithm>
#include <iostream>
#include <cassert>
#include <thread>
#include <atomic>
#include <functional>

class MultithreadedAligner
{
public:
  static const int ctThreads = 8;

  struct ThreadData {
    struct RowData {
      int row;
      int valBestInCol;
      char y;
    };

    ThreadData() {}
    ThreadData(ThreadData &td) {}
    vector<int> _col;
    int *col;
    vector<RowData> localRowData;
    atomic<int> iRow;
    EndPoint best;
    function<void(void)> f = nullptr;
    atomic<int> go = false;
    atomic<int> die = false;
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

    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      tds[iThread].f = [iThread, rev, iRow0, iRow1, iCol0, iCol1,
        rowStart, row, &valBestInRow, &valBestInCol, &x, &y, &sc, this]()
      {
        const Scoring scLocal = sc;
        EndPoint best;

        const int dir = rev ? -1 : 1;

        int iThCol0 = iCol0 + (iCol1 - iCol0) * (iThread - 1) / ctThreads;
        int iThCol1 = iCol0 + (iCol1 - iCol0) * iThread / ctThreads - (iThread == ctThreads ? 0 : dir);

        int iThRow0 = iRow0;
        int iThRow1 = iRow1 + dir;
        iThCol1 += dir;

        const int cpyStart = min(iThCol0, iThCol1 - dir);
        const int cpyCount = abs(iThCol1 - iThCol0);

        auto *localRowData = tds[iThread].localRowData.data();
        for (int iCol = cpyStart; iCol < cpyStart + cpyCount; ++iCol) {
          localRowData[iCol].row = row[iCol];
          localRowData[iCol].valBestInCol = valBestInCol[iCol];
          localRowData[iCol].y = y[iCol + rev];
        }

        tds[iThread].col[iThRow0 - dir] = localRowData[iThCol1 - dir].row;
        for (int iRow = iThRow0; iRow != iThRow1; iRow += dir) {
          while (dir * tds[iThread - 1].iRow < dir * iRow) {
            this_thread::yield();
          }
          int prevColRow = tds[iThread - 1].col[iRow - dir];
          int prevCol = tds[iThread - 1].col[iRow];

          const char xLocal = x[iRow + rev];
          int valBestInRowLocal = valBestInRow[iRow];
          for (int iCol = iThCol0; iCol != iThCol1; iCol += dir) {
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
          tds[iThread].col[iRow] = localRowData[iThCol1 - dir].row;
          tds[iThread].iRow = iRow;
        }
        for (int iCol = cpyStart; iCol < cpyStart + cpyCount; ++iCol) {
          row[iCol] = localRowData[iCol].row;
          valBestInCol[iCol] = localRowData[iCol].valBestInCol;
        }
        tds[iThread].best.Add(best);
      };
      tds[iThread].go = 1;
    }
    EndPoint r;
    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      while (tds[iThread].go) {
        this_thread::yield();
      }
      r.Add(tds[iThread].best);
    }
    return r;
  }

  void init(int ctRow, int ctCol)
  {
    tds.resize(ctThreads + 1);
    for (auto &td : tds) {
      td._col.resize(ctRow + 3);
      td.col = td._col.data() + 1;
      td.localRowData.resize(ctCol + 2);
    }
    ts.resize(ctThreads + 1);
    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
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
    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      tds[iThread].die = 1;
    }
    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      ts[iThread].join();
    }
  }
};
