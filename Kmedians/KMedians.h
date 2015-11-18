#pragma once

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <vector>
#include <iostream>
#include <list>
#include <map>
#include <algorithm>    // std::random_shuffle
#include <ctime>        // std::time
#include <cstdlib>      // std::rand, std::srand
#include <set>

//______________________________________________________________________________
//complementary functions declaration...........................................
int		kmed_hammingDistanceInt(const int x, const int y);
//______________________________________________________________________________
// ORIGINAL CLASS, DECALRATION VARIABLE......................................
/*private static final long serialVersionUID = -253576488350574703L;

	private int nbThread = 1;
	
	int[][] points;
	double[][] centers;
	double[] meanDistance;
	int[] populationInCluster;
	int[] center_Point;
	ArrayList<ArrayList<Integer>> pointsInCluster;
	
	int nbCenters;
	int dimension = -1;
	
	int pointAssignedToCenter[];
	double distanceMatrix[][];
	boolean hasMoved[];
	boolean nullWord;
	
	int maxIterations = 10000;
	
	// random generator
	Random ran;
*/

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
	cv::RNG				ran;
	//Constructor destructor....................................................
	KMedians(){}
	~KMedians(){}
	//Main funtions.............................................................
	void				init_centers();
	int					makeAssigment();
	void				computeDistance();

};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////main functions//////////////////////////////////////////////


/*
	private void initCenters()
	{
		pointsInCluster = new ArrayList<ArrayList<Integer>>(nbCenters);
		dimension = points[0].length;
		centers = new double[nbCenters][dimension];
		meanDistance = new double[nbCenters];
		populationInCluster = new int[centers.length];
		//center_Point = new int[centers.length];
		//Arrays.fill(center_Point, -1);
		Arrays.fill(populationInCluster, 0);
		
		//random assignment
		System.out.println("Initializing centers...");
		pointAssignedToCenter = new int[points.length];
		
		//-1 == no assignment
		Arrays.fill(pointAssignedToCenter, -1);
		List<Integer> listRandom = new ArrayList<Integer>(); 
		
		//pick a random point for each cluster
		for(int i = 0 ; i < centers.length; i++)
		{
			int indexPoint = ran.nextInt(points.length);
			while (listRandom.contains(indexPoint))
				indexPoint = ran.nextInt(points.length);
			listRandom.add(indexPoint);
			pointAssignedToCenter[indexPoint] = i;
			
			/* // --> original commented code
			ArrayList<Integer> p =  new ArrayList<Integer>();
			p.add(indexPoint);
			pointsInCluster.set(i, p);
			*/ // --> original commented code
			
			/*populationInCluster[i]++;
			//center_Point[i] = indexPoint;
			if(i%(centers.length/20 + 1) == 0)
				System.out.print(".");
		}
		System.out.println();
		
		//distance matrix and has moved
		distanceMatrix = new double[points.length][centers.length];
		hasMoved = new boolean[centers.length];
		Arrays.fill(hasMoved, true);
		
		System.out.println("Centers randomly initialized.");
	}
*/

void KMedians::init_centers(){
	pointsInCluster_.resize(nbCenters_);
	dimension_	= points_.cols;
	centers_.create(nbCenters_, dimension_);
	meanDistance_.resize(nbCenters_);
	populationInCluster_.resize(centers_.rows, 0);
	//_________________________________________________________________________
	//random assiggmnet 
	std::cout << "Initializing centers...\n";
	pointAssignedToCenter_.resize(points_.rows,-1);			//-1 = no assig.
	std::set<int>	listRandom;
	//_________________________________________________________________________
	//pick a random point for each cluster
	for (int i = 0; i < centers_.rows; ++i)
	{
		int indexPoint = ran.uniform(0, points_.rows);
		while ( listRandom.find(indexPoint) != listRandom.end() )
			indexPoint = ran.uniform(0, points_.rows);
		listRandom.insert(indexPoint);
		pointAssignedToCenter_[indexPoint] = i;
		++populationInCluster_[i];
		if (i % (centers_.rows / 20 + 1) == 0)
            std::cout << "." << std::endl;
	}
	distanceMatrix_.create(points_.rows, centers_.rows);
	hasMoved_.resize(centers_.rows, true);
	std::cout << "\nCenters randomly initialized\n";
}

////////////////////////////////////////////////////////////////////////////////

/*private int makeAssignment() {

		int nbm = 0;
		
		Arrays.fill(populationInCluster, 0);
		Arrays.fill(hasMoved, false);
		
		//pointsInCluster = new ArrayList<ArrayList<Integer>>(nbCenters);
		
		for(int i = 0 ; i < points.length; i++)
		{
			//find distance min
			int indexMin = 0;
			double dist = distanceMatrix[i][0];
			for(int m = 0 ; m < centers.length; m++)
			{
				if(distanceMatrix[i][m] < dist)
				{
					dist = distanceMatrix[i][m];
					indexMin = m;
				}
			}
			
			populationInCluster[indexMin]++;
			
			//compare to original
			int oldIndex = pointAssignedToCenter[i];
			if(oldIndex != indexMin)
			{
				//got one more move
				nbm++;
				//centers will change
				hasMoved[indexMin] = true;
				if(oldIndex != -1)
					hasMoved[oldIndex] = true;
				
				//make assignment
				pointAssignedToCenter[i] = indexMin;
				
				/*=========original comment
				ArrayList<Integer> a = pointsInCluster.get(indexMin);
				a.add(i);
				pointsInCluster.set(indexMin, a);
				*/ //=========original comment
				
			/*}
		}	
		return nbm;
	}
	*/

int KMedians::makeAssigment()
{
	int nbm = 0;
	std::fill(populationInCluster_.begin(), populationInCluster_.end(), 0);
	std::fill(hasMoved_.begin(), hasMoved_.end(), false);
	for (int i = 0; i < points_.rows; ++i)
	{
		double dist		= distanceMatrix_(i, 0);
		int	indexmin	= 0;
		for (int m = 0; m < centers_.rows; ++m)
		{
			if (distanceMatrix_(i, m) < dist)
			{
				dist	 = distanceMatrix_(i, m);
				indexmin = m;
			}
		}
		//compare with the original
		++populationInCluster_[indexmin];
		int oldIndex = pointAssignedToCenter_[i];
		if (oldIndex != indexmin)
		{
			//got one more move
			++nbm;
			//centers will change
			hasMoved_[indexmin] = true;
			if (oldIndex != -1)
				hasMoved_[oldIndex] = true;
			//make assigment
			pointAssignedToCenter_[i] = indexmin;
		}
		
	}
	return nbm;
}

////////////////////////////////////////////////////////////////////////////////
/*private void computeDistance() {
		
		ExecutorService es = Executors.newFixedThreadPool(nbThread);

		List<Callable<Object>> tasks = new ArrayList<Callable<Object>>();
		
		for(int n = 0 ; n < points.length; n++)
		{
			final int i = n;

			Callable<Object> r = new Callable<Object>(){

				public Object call()
				{

						for(int c = 0 ; c < centers.length; c++)
						{
							//if center didn't move, continue
							if(!hasMoved[c])
								continue;
							
							//if hasn't points in it, continue;
							if(populationInCluster[c] == 0)
							{
								distanceMatrix[i][c] = Double.POSITIVE_INFINITY;
								continue;
							}

							double distance = 0;
							int[] p = points[i];
							double[] center = centers[c];
							
							for(int d = 0; d < dimension; d++) {
								
								//Alterado: Verifica se o centro � NaN pois, se n�o verificar, na tranforma��o para int o NaN vira Zero e afeta no pr�ximo if
								if(Double.isNaN(center[d])) {
									distance = Double.NaN;
									break; //break pois se o primeiro valor do centro ja � NaN o resto tamb�m � NaN
								}
								else {
									distance += hammingDistance( (p[d]& 0xff), (((int)center[d]) & 0xff) );	
								}
							}
							    
							if(Double.isNaN(distance)) //afeta neste if
								distance = Double.POSITIVE_INFINITY;
							distanceMatrix[i][c] = distance; 
						}
						return null;
				};
			};
			
			tasks.add(r);
		}
		
		try {
			//waiting for all threads to finish
			es.invokeAll(tasks);
			
		} catch (InterruptedException e) {
			System.err.println("Fatal Error : Main thread interreupted!");
			e.printStackTrace();
			System.exit(-1);
		}
	}
*/

void KMedians::computeDistance()
{

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
