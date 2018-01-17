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

    ThreadData() {}
    ThreadData(const ThreadData &td) {}
    vector<int> _col;
    int *col;
    vector<RowData> localRowData;
    atomic<int> iRow;
    EndPoint best;
    function<void(void)> f = nullptr;
    atomic<int> go{false};
    atomic<int> die{false}; 

  };

  vector<ThreadData> tds;

  vector<thread> ts;

  EndPoint _alignMultithreaded(int ctRow, int ctCol, const char *x, const char *y, int rowStart,
      int *row, vector<int> &valBestInRow, vector<int> &valBestInCol, const Scoring &sc)
  {
    const int ctThreads = min(ctMaxThreads, ctCol);

    tds[0].col[-1] = row[-1];
    for (int iRow = 0; iRow < ctRow; ++iRow) {
      tds[0].col[iRow] = rowStart;
    }
    tds[0].iRow = ctRow - 1;

    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      tds[iThread].iRow = -1;
      tds[iThread].best = EndPoint{};
    }

    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      tds[iThread].f = [iThread, ctThreads, ctRow, ctCol,
        rowStart, row, &valBestInRow, &valBestInCol, x, y, &sc, this]()
      {
        const Scoring scLocal = sc;
        EndPoint best;
        int iThCol0 = ctCol * (iThread - 1) / ctThreads;
        int iThCol1 = ctCol * iThread / ctThreads;

        auto *localRowData = tds[iThread].localRowData.data();
        for (int iCol = iThCol0; iCol < iThCol1; ++iCol) {
          localRowData[iCol].row = row[iCol];
          localRowData[iCol].valBestInCol = valBestInCol[iCol];
          localRowData[iCol].y = y[iCol];
        }
        tds[iThread].col[-1] = row[iThCol1 - 1];

        for (int iRow = 0; iRow < ctRow; ++iRow) {
          while (tds[iThread - 1].iRow < iRow) {
            this_thread::yield();
          }
          int prevColRow = tds[iThread - 1].col[iRow - 1];
          int prevCol = tds[iThread - 1].col[iRow];

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
        tds[iThread].best.Add(best);
      };
      tds[iThread].go = 1;
    }
    EndPoint r;
    for (int iThread = 1; iThread <= ctThreads; ++iThread) {
      while (tds[iThread].go) {
        this_thread::yield();
      }
      const EndPoint &b = tds[iThread].best;
      if (b.val > r.val || b.val == r.val &&
          (b.p.x < r.p.x || b.p.x == r.p.x && b.p.y < r.p.y)) {
          r = b;
      }
    }
    return r;
  }

  void init(int ctRow, int ctCol)
  {
    ctMaxThreads = min(ctCol, 16);

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
  }
};
