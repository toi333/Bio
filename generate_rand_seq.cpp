#include <iostream>
#include <stdlib.h>
#include <string>

using namespace std;

int main(int argc, char**argv) {
  if (argc == 1) {
    cout << "Give size of sequance" << endl;
    return 1;
  }

  srand(time(0));

  long seq_len = atol(argv[1]);

  while(seq_len > 0) {
    seq_len -= 80;
    cout << ">" << seq_len << endl;
    for(int i = 0; i < 80; ++i) {
      cout <<  "ACTG"[rand() % 4];
    }
    cout << endl;
  }


  return 0;
}
