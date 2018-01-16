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

void test(const string &s1, const string &s2)
{
    HirschbergAA algo;
    //SquareAA algo;
    //NaiveCubeAA algo;
    Scoring sc;

    int start = clock();
    Alignment sol = algo.align(s1, s2, sc);
    float runTimeAlg = (float)(clock() - start) / CLOCKS_PER_SEC;
    const int validationScore = sol.calcScore(s1, s2, sc);
    if (validationScore != sol.score) {
        cout << "VALIDATION FAILED!" << endl;
        cout << "Score should be: " << validationScore << endl;
        sol.output(cout, s1, s2);
    }
    cout << "Runtime algo:   " << runTimeAlg << endl;

    //SquareAA algoSq;
    ////int startSq = clock();
    //Alignment sol2 = algoSq.align(s1, s2, sc);
    ////float runTimeSq = (float)(clock() - startSq) / CLOCKS_PER_SEC;

    ////cout << "Runtime square: " << runTimeSq << endl;
    //if (sol.score != sol2.score) {
    //    cout << "FAILED!" << endl;
    //    cout << "H: " << sol.score << endl;
    //    sol2.output(cout, s1, s2);
    //} else {
    //    cout << "OK!" << endl;
    //}
    cout << endl;
}

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

  int cutoff = 30000;
  for (auto &s : vs) {
    s = s.substr(0, cutoff);
  }

  cout << "TEST: " << file << ' ' << i1 << ' ' << i2 << endl;

  test(vs[i1], vs[i2]);
}

int main()
{
  srand(time(0));

  init();

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
  //test("res/streptococcus_references.fasta", 0, 1);

  for (int l = 10; l < 1000; ++l) {
      for (int t = 0; t < 10; ++t) {
          string s1, s2;
          for (int i = 0; i < l; ++i) {
              s1.push_back("ACTG"[rand() % 4]);
              s2.push_back("ACTG"[rand() % 4]);
          }
          test(s1, s2);
      }
  }

  system("pause");
}