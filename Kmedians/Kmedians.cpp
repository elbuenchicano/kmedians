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

	vector<int> a;
	a.resize(3, 5);
	copy ( a.begin(), a.end(), ostream_iterator<int> (std::cout," "));

	RNG rng( 0xFFFFFFFF ), rng2;

	rng2 = rng;

	return 0;
}
