#include <math.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>

int main(){
  //std::string filename;
  char filename[32];
  int i=4;
  sprintf(filename, "output%03d.txt", i);
  std::cout<<filename<<std::endl;
  //filename = "Dense"+std::to_string(4)t+".txt";
  //std::cout<<filename<<std::endl;
}
