#pragma once

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <iostream>
#include <list>
#include <map>

//______________________________________________________________________________
//complementary functions declaration...........................................
int kmed_hammingDistanceInt(const int x, const int y);
//______________________________________________________________________________

struct KMedians
{
	int					nbCenters_,
						dimension_		= -1,
						maxIterations_	= 10000;
	
	bool				nullWord_;

	std::vector<bool>	hasMoved_;
	cv::Mat_<int>		points_;
	cv::Mat_<float>		centers_,
						distanceMatrix_;
	
	std::vector<double>	meanDistance_;
	std::vector<int>	populationInCluster_,
						centerPoint_,
						pointAssignedToCenter_;
	using ArrayListInt  = std::vector<std::vector<int> >;
	ArrayListInt		pointsInCluster_;
	cv::RNG				rng;
	//main functions............................................................
	KMedians(){}
	~KMedians(){}
	void				init_centers();

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////main functions//////////////////////////////////////////////

void KMedians::init_centers(){
	pointsInCluster_.resize(nbCenters_);
	dimension_	= points_.cols;
	centers_.create(nbCenters_, dimension_);
	meanDistance_.resize(nbCenters_);
	populationInCluster_.resize(centers_.rows, 0);
	//_________________________________________________________________________
	//random assiggmnet 
	std::cout << "Initializing centers..." << std::endl;
	pointAssignedToCenter_.resize(points_.rows,-1);			//-1 = no assig.
	std::vector<int>	listRandom;
	//_________________________________________________________________________
	//pick a random point for each cluster
	for (int i = 0; i < centers_.rows; ++i){
		int indexPoint = rng.uniform(0, centers_.rows);

	}
}



////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////Complementary functions/////////////////////////////////////
//simple binary hamming distance between two integers
int kmed_hammingDistanceInt(const int x, const int y){
	int dist = 0;
	for ( int val = (x ^ y); val; ++dist ) val &= val - 1;		
	return dist;
}
