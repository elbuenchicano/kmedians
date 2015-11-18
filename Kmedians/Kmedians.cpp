// Kmedians.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include "KMedians.h"
#include <iterator>


using namespace cv;
using namespace std;



int main(int argc, char** argv)
{
	/*Mat image;
	image = imread("d:/chela.jpg", CV_LOAD_IMAGE_COLOR);   // Read the file

	if (!image.data)                              // Check for invalid input
	{
	cout << "Could not open or find the image" << std::endl;
	return -1;
	}

	namedWindow("Display window", WINDOW_AUTOSIZE);// Create a window for display.
	imshow("Display window", image);                   // Show our image inside it.

	waitKey(0);                                          // Wait for a keystroke in the window*/

	/*vector<int> a;
	a.resize(3, 5);
	copy ( a.begin(), a.end(), ostream_iterator<int> (std::cout," "));

	RNG rng( 0xFFFFFFFF ), rng2;

	rng2 = rng;
	cout<< rng2.uniform(0,1000);/**/

	/*std::srand ( unsigned ( std::time(0) ) );
	std::vector<int> myvector;

	// set some values:
	for (int i=1; i<10; ++i) myvector.push_back(i); // 1 2 3 4 5 6 7 8 9

	// using built-in random generator:
	std::random_shuffle ( myvector.begin(), myvector.end() );

	// using myrandom:
	std::random_shuffle ( myvector.begin(), myvector.end(), myrandom);

	// print out content:
	std::cout << "myvector contains:";
	for (std::vector<int>::iterator it=myvector.begin(); it!=myvector.end(); ++it)
	std::cout << ' ' << *it;

	std::cout << '\n';
	*/

	set<int> st{ 2, 6, 2, 7, 1 };
	copy ( st.begin(), st.end(), ostream_iterator<int> (std::cout," "));
	//auto it = st.find(2);
	if (st.find(3)!=st.end())cout << "ok";

	return 0;
}
