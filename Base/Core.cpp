#include <Core.h>

#include <iostream>
#include <algorithm>
#include <cstring>

char nidToIdx[300];

void init()
{
  memset(nidToIdx, -1, sizeof nidToIdx);
  nidToIdx['A'] = 0;
  nidToIdx['C'] = 1;
  nidToIdx['G'] = 2;
  nidToIdx['T'] = 3;
}


vector<char> encode(const string &a)
{
  vector<char> x;
  x.resize(a.size() + 2);
  x[0] = 4;
  for (int i = 0; i < (int)a.size(); ++i) {
    x[i + 1] = nidToIdx[a[i]];
  }
  x[a.size() + 1] = 4;
  return x;
}

void Alignment::output(ostream &os, const string &a, const string &b) const
{
  os << score << endl;

  if (matches.empty()) {
    os << "Null alignment" << endl;
    return;
  }

  string sx, sy;
  Vec2i cur = matches[0];
  os << cur.x << ' ' << cur.y << endl;
  for (int i = 1; i < (int)matches.size(); ++i) {
    const Vec2i next = matches[i];
    if (next.x == cur.x + 1 && next.y == cur.y + 1) {
      cur = next;
      sx += a[cur.x - 1];
      sy += b[cur.y - 1];
      if (sx.size() == 80) {
        os << sx << endl;
        os << sy << endl;
        os << endl;
        sx = "";
        sy = "";
      }
    } else {
      while (cur.x < next.x) {
        ++cur.x;
        sx += a[cur.x - 1];
        sy += '-';
        if (sx.size() == 80) {
          os << sx << endl;
          os << sy << endl;
          os << endl;
          sx = "";
          sy = "";
        }
      }
      while (cur.y < next.y) {
        ++cur.y;
        sx += '-';
        sy += b[cur.y - 1];
        if (sx.size() == 80) {
          os << sx << endl;
          os << sy << endl;
          os << endl;
          sx = "";
          sy = "";
        }
      }
    }
  }
  if (sx.size() > 0) {
    os << sx << endl;
    os << sy << endl;
    os << endl;
  }
}

int Alignment::calcScore(const string &a, const string &b, const Scoring &sc) const
{
  int validationScore = 0;

  vector<char> x = encode(a);
  vector<char> y = encode(b);

  Vec2i prev = matches[0];
  for (int i = 1; i < (int)matches.size(); ++i) {
    const Vec2i cur = matches[i];
    if (cur.x == prev.x && cur.y == prev.y) {
      cerr << "FAIL" << endl;
      return -1;
    }
    if (cur.x < prev.x || cur.y < prev.y) {
      cerr << "FAIL" << endl;
      return -1;
    }
    if (cur.x != prev.x && cur.y != prev.y) {
      if (cur.x != prev.x + 1 || cur.y != prev.y + 1) {
        cerr << "FAIL" << endl;
        return -1;
      }
      validationScore += sc.match(x[cur.x], y[cur.y]);
    } else {
      validationScore -= sc.b + sc.k * (max(cur.x - prev.x, cur.y - prev.y) - 1);
    }
    prev = cur;
  }

  return validationScore;
}

void Alignment::compress()
{
  if (matches.size() <= 2) {
    return;
  }

  vector<Vec2i> old;
  old.swap(matches);
  matches.push_back(old[0]);
  for (int i = 1; i + 1 < old.size(); ++i) {
    bool diff =
      ((old[i].x != old[i - 1].x) || (old[i].x != old[i + 1].x)) &&
      ((old[i].y != old[i - 1].y) || (old[i].y != old[i + 1].y));
    if (diff) {
      matches.push_back(old[i]);
    }
  }
  matches.push_back(old.back());
}
