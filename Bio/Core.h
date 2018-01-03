#pragma once

#include <fstream>
#include <vector>
#include <string>

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

typedef Mat2<int> Mat2i;

class Scoring
{
public:
  static const int mp = 3;
  //static const int mn = -4;
  static const int mn = -100;
  static const int ni = -(1 << 28);

  int match(int a, int b) const {
    return a == b ? mp : mn;
  }

  int m[5][5] = {
    { mp, mn, mn, mn, ni },
    { mn, mp, mn, mn, ni },
    { mn, mn, mp, mn, ni },
    { mn, mn, mn, mp, ni },
    { ni, ni, ni, ni, ni },
  };

  int k = 1;
  int b = 4;
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
