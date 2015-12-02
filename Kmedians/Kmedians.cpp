// Kmedians.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "KMedians.h"
#include <iterator>
#include <thread>
#include <string.h>
#include <fstream>
#include <sstream>

using namespace cv;
using namespace std;



void loadTest(string file, Mat_<int> & des){
  fstream arc(file);
  string  line;
  char    buffer[20];
  if (!arc.is_open()) return;
  while (!arc.eof()){
    getline(arc, line);
    auto pos  = line.find(">");
    line = line.substr(pos+2);
    Mat_<int> row(1, 16);
    row = row * 0;
    for (int i = 0, posLine = 0; posLine < line.size() && i < 16; ++posLine){
      if (line[posLine] == ' ') ++i;
      else{
        row(0, i) = row(0, i) * 10;
        row(0, i) = row(0, i) + line[posLine] - 48;
      }
    }
    des.push_back(row);
  }
  arc.close();
}


int main(int argc, char** argv)
{
  Mat_<int> des;
  loadTest("d:/kmedinput.txt", des);
  RNG       ran( 0 );
  KMedians  kmed(des, 10, 5000, 4, ran);
  FileStorage fs("d:/outcenters.yml", FileStorage::WRITE);
  fs << "centers" << kmed.centers_;
	return 0;

}
