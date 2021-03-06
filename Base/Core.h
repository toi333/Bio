#pragma once

#include <fstream>
#include <vector>
#include <string>
#include <climits>

using namespace std;

extern char nidToIdx[300];

void init();

class Vec2i {
  public:
    int x, y;
};

template<class T>
struct Mat2 {
  public:
    int n, m;
    vector<T> t;

    void resize(int n, int m) {
      this->n = n;
      this->m = m;
      t.resize(n * m);
    }

    T operator[](const Vec2i &v) const {
      return t[m * v.x + v.y];
    }

    T& operator[](const Vec2i &v) {
      return t[m * v.x + v.y];
    }
};

struct EndPoint
{
  int val = INT_MIN;
  Vec2i p = Vec2i{ -1, -1 };

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

typedef Mat2<int> Mat2i;

class Scoring
{
  public:
    const int mp = 1;
    // const int _match = 1;
    // const int mp = 3;
    //const int mn = -4;
    const int mn = -3;
    // const int mismatch = -3;
    const int ni = -(1 << 28);

    int match(int a, int b) const {
      return a == b ? mp : mn;
    }

    // const int gapOpen = 5;
    const int b = 5;
    // const int gapExtend = 2;
    const int k = 2;
};

class Alignment {
  public:
    int score = 0;
    vector<Vec2i> matches;

    void output(ostream &is, const string &a, const string &b) const;

    int calcScore(const string &a, const string &b, const Scoring &sc) const;

    void compress();
};

class AlignmentAlgorithm {
  public:
    virtual ~AlignmentAlgorithm() {}

    virtual Alignment align(const string &a, const string &b, const Scoring &sc) = 0;
};

vector<char> encode(const string &a);
