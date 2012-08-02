#include "ttfmm.hpp"

using namespace std;
int main(int argc, char* argv[]){
  char *filename = argv[1];
  double alpha = atof(argv[2]);
  double lambda = atof(argv[3]);
  int k_max = atoi(argv[4]);
  
  MM m(alpha, lambda, k_max);
  ifstream ifs;
  ifs.open(filename, ios::in);
  string line;
  while(getline(ifs, line)){
    // file format
    // timestamp \t w_1 \t w_2 ...
    vector<string> elem = split_string(line, "\t");
    m.set_document(elem);
  }
  ifs.close();

  m.fitting();
}
