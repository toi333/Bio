#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <string>
#include <ctime>

#include <Base/Core.h>
#include <Bio/Algorithms/NaiveCubeAlgo.h>
#include <Bio/Algorithms/HirschbergAlgo.h>
#include <Bio/Algorithms/SquareAlgo.h>

using namespace std;

const vector<string> readFastaVector(const char *file)
{
  ifstream fi(file);
  vector<string> vs;
  string l;
  while (getline(fi, l)) {
    if (l[0] == '>') {
      vs.push_back(string());
    } else {
      vs.back() += l;
    }
  }
  fi.close();

  return vs;
}

const string readFastaString(const char *file)
{
  auto vs = readFastaVector(file);
  string ret = "";
  long totalSize = 0;

  for(const auto& s: vs) {
    totalSize += s.size();
  }
  ret.reserve(totalSize);

  for(const auto& s: vs) {
    ret += s;
  }

  return ret;
}


