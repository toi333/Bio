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

extern int testcuda();

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

  int cutoff = 50000;
  for (auto &s : vs) {
    s = s.substr(0, cutoff);
  }

  HirschbergAA algo;
  //SquareAA algo;
  //NaiveCubeAA algo;
  Scoring sc;

  cout << "TEST: " << file << ' ' << i1 << ' ' << i2 << endl;

  int start = clock();
  Alignment sol = algo.align(vs[i1], vs[i2], sc);
  float runTimeAlg = (float)(clock() - start) / CLOCKS_PER_SEC;
  const int validationScore = sol.calcScore(vs[i1], vs[i2], sc);
  if (validationScore != sol.score) { 
    cout << "VALIDATION FAILED!" << endl;
    cout << "Score should be: " << validationScore << endl;
    sol.output(cout, vs[i1], vs[i2]);
  }
  cout << "Runtime algo:   " << runTimeAlg << endl;

  //SquareAA algoSq;
  //int startSq = clock();
  //Alignment sol2 = algoSq.align(vs[i1], vs[i2], sc);
  //float runTimeSq = (float)(clock() - startSq) / CLOCKS_PER_SEC;

  //cout << "Runtime square: " << runTimeSq << endl;
  //if (sol.score != sol2.score) {
  //  cout << "FAILED!" << endl;
  //  cout << "H: " << sol.score << endl;
  //  sol2.output(cout, vs[i1], vs[i2]);
  //} else {
  //  cout << "OK!" << endl;
  //}
  cout << endl;
}

int main()
{
  init();

  testcuda();

  //test("res/test.fasta", 0, 0);
  //test("res/test.fasta", 0, 1);
  //test("res/test.fasta", 0, 2);
  //test("res/test.fasta", 3, 4);
  //test("res/test.fasta", 5, 6);
  //test("res/test.fasta", 6, 5);
  //test("res/test.fasta", 7, 8);
  //test("res/test.fasta", 9, 10);
  //test("res/test.fasta", 9, 11);
  //test("res/test.fasta", 12, 13);
  //test("res/test.fasta", 14, 15);
  //test("res/test.fasta", 15, 14);
  //test("res/test.fasta", 16, 17);
  //test("res/test.fasta", 18, 18);
  //test("res/test.fasta", 19, 20);
  //test("res/test.fasta", 20, 19);
  //test("res/test.fasta", 21, 22);
  //test("res/test.fasta", 23, 24);
  test("res/streptococcus_references.fasta", 0, 1);

  system("pause");
}