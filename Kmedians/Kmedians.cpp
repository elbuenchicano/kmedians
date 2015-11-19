// Kmedians.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "KMedians.h"
#include <iterator>
#include <thread>

using namespace cv;
using namespace std;

//std::cout << "Hello from thread " << std::this_thread::get_id() << std::endl;
struct X
{

	vector<int> mX;

	void F(int i, int f)
	{
		for (int it = i; it < f; ++it)
			cout << mX[it] << " ";
	}
	X()
	{
		mX.resize(20);
		for (int i = 0; i < 20; ++i)
			mX[i] = i;
	}

	void G()
	{
		vector<thread> vt;
		for (int i = 0; i < 20; i += 4){
			vt.push_back(thread ( &X::F , this, i, i + 4) );
		}

		for (auto & it: vt){
			it.join();
		}
	}

};
 
void thread_sample()
{
	X x;
	x.G();
}

int main(int argc, char** argv)
{
	//thread_sample();
	Mat_<int> x(1,1);
	x(0, 0) = NAN;
	cout << x(0,0);
	return 0;
}
