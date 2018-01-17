#include <iostream>
#include <string>
#include <ctime>
#include <Bio/Util.hpp>


using namespace std;

int main(int argc, char**argv) 
{
  init();
  if(argc < 3) {
    cout << "No input files" << endl;
    cout << "./exe file1.fasta file2.fasta" << endl;
    
    return 1;
  }

  const string genom1= readFastaString(argv[1]);
  const string genom2= readFastaString(argv[2]);
  HirschbergAA algo;
  Scoring sc;
  int start = clock();
  Alignment sol = algo.align(genom1, genom2, sc);
  float runTimeAlg = (float)(clock() - start) / CLOCKS_PER_SEC;
  cout << "Algo runtime: " << runTimeAlg << endl;
  cout << "Score: " << sol.score << endl;


# ifndef __linux__
  system("pause");
#endif

  return 0;
}
