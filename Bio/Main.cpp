#include <algorithm>
#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <string>
#include <ctime>

#include <Core.h>
#include <Algorithms/NaiveCubeAlgo.h>
#include <Algorithms/HirschbergAlgo.h>
#include <Algorithms/SquareAlgo.h>

using namespace std;

void test(const char *file, int i1, int i2)
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

  int cutoff = 1000;
  for (auto &s : vs) {
    s = s.substr(0, cutoff);
  }

  HirschbergAA algo;
  //SquareAA algo;
  //NaiveCubeAA algo;
  Scoring sc;

  int start = clock();
  Alignment sol = algo.align(vs[i1], vs[i2], sc);
  float runTimeAlg = (float)(clock() - start) / CLOCKS_PER_SEC;
  //cout << "Run time: " << (float)(clock() - start) / CLOCKS_PER_SEC << endl;
  const int validationScore = sol.calcScore(vs[i1], vs[i2], sc);
  if (validationScore != sol.score) {
    cout << "VALIDATION FAILED: " << file << ' ' << i1 << ' ' << i2 << endl;
    cout << "Score should be: " << validationScore << endl;
    sol.output(cout, vs[i1], vs[i2]);
  }

  //sol.output(cout, vs[i1], vs[i2]);

  SquareAA algoSq;
  int startSq = clock();
  Alignment sol2 = algoSq.align(vs[i1], vs[i2], sc);
  float runTimeSq = (float)(clock() - startSq) / CLOCKS_PER_SEC;
  //const int validationScore = sol2.calcScore(vs[i1], vs[i2], sc);
  //if (validationScore != sol2.score) {
  //  cout << "VALIDATION FAILED: " << file << ' ' << i1 << ' ' << i2 << endl;
  //  cout << "Score should be: " << validationScore << endl;
  //  sol2.output(cout, vs[i1], vs[i2]);
  //}

  if (sol.score != sol2.score) {
    cout << "FAILED: " << file << ' ' << i1 << ' ' << i2 << endl;
    cout << "H: " << sol.score << endl;
    sol2.output(cout, vs[i1], vs[i2]);
  } else {
    cout << "OK: " << file << ' ' << i1 << ' ' << i2 << endl;
  }
  cout << "Runtime algo:   " << runTimeAlg << endl;
  cout << "Runtime square: " << runTimeSq << endl;
  cout << endl;
}

int main()
{
  init();

  test("res/test.fasta", 0, 0);
  test("res/test.fasta", 0, 1);
  test("res/test.fasta", 0, 2);
  test("res/test.fasta", 3, 4);
  test("res/test.fasta", 5, 6);
  test("res/test.fasta", 6, 5);
  test("res/test.fasta", 7, 8);
  test("res/test.fasta", 9, 10);
  test("res/test.fasta", 9, 11);
  test("res/test.fasta", 12, 13);
  test("res/test.fasta", 14, 15);
  test("res/test.fasta", 15, 14);
  test("res/test.fasta", 16, 17);
  test("res/test.fasta", 18, 18);
  test("res/streptococcus_references.fasta", 0, 1);

  system("pause");
}