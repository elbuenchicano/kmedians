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
#include <mutex>
#include <thread>

//______________________________________________________________________________
//complementary functions declaration...........................................
int		kmed_hammingDistanceInt(const int x, const int y);
template<class t> void	fillMatRow(cv::Mat_<t> & src, int row, t val);
//______________________________________________________________________________
// ORIGINAL CLASS, DECALRATION VARIABLE......................................
/*private static final long serialVersionUID = -253576488350574703L;

	private int nbThread_ = 1;
	
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
	int       					nbCenters_,
						          dimension_		  = -1,
						          maxIterations_	= 10000,
						          nbThread_		    = 1;
	
	bool				        nullWord_;

	std::vector<bool>	  hasMoved_;
	cv::Mat_<int>		    points_;
	cv::Mat_<float>		  centers_,
						          distanceMatrix_;
	
	std::vector<double>	meanDistance_;
	std::vector<int>	  populationInCluster_,
						          centerPoint_,
						          pointAssignedToCenter_;
	using ArrayListInt  = std::vector<std::vector<int> >;
	ArrayListInt		    pointsInCluster_;
	cv::RNG				      ran_;
	//Constructor destructor....................................................
  KMedians  ( cv::Mat_<int>, int, int, int, cv::RNG );
	KMedians  (){}
	~KMedians (){}
	//Main funtions.............................................................
	void				init_centers	();
	int					makeAssigment	();
	void				computeDistance	();
	void				computeCenters	();
	void				computeCentersMedian	();
	//Support functions.........................................................
	void				computeDistanceFor(int, int);
  
};
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////main functions//////////////////////////////////////////////

KMedians::KMedians(cv::Mat_<int> points, int nbCenters, int maxIterations, int nbThreads, cv::RNG r) :
points_(points),
nbCenters_(nbCenters),
maxIterations_(maxIterations),
nbThread_(nbThreads),
ran_(r)
{

  init_centers();
  std::cout << "Starting clustering with " << nbThread_ << "threads.\n";
  int nbMoves   = 0,
      iteration = 1;
  do{
    computeCentersMedian();
    //compute distance
    computeDistance();
    //number of moving centers
    int nbmc = 0;
    for (int c = 0; c < centers_.rows; ++c){
      if (hasMoved_[c])
        ++nbmc;
    }
    std::cout << "iteration " << nbMoves << " points moved " << nbmc << " centers moved.\n";
  } while (nbMoves != 0 && ++iteration < maxIterations_);
  std::cout << "\n Clustering done, cleaning ....\n";
  std::fill(meanDistance_.begin(), meanDistance_.end(), 0);
  for (int c = 0; c < centers_.rows; ++c){
    int nbp = 0;
    for (int i = 0; i < points_.rows; ++i){
      if (pointAssignedToCenter_[i] == c){
        meanDistance_[c] += distanceMatrix_(i,c);
        ++nbp;
      }
    }
    if (nbp > 0)
      meanDistance_[c] /= (double)nbp;
  }
  //cleaning empty clusters 
  std::vector< cv::Mat_<float> >  listOfCenters;
  std::vector< double >       listOfMeanDist;
  std::vector< int >          listOfPopulation;
  for (int c = 0; c < centers_.rows; ++c){
    listOfCenters.push_back(centers_.row(c));
    listOfMeanDist.push_back(meanDistance_[c]);
    listOfPopulation.push_back(populationInCluster_[c]);
  }
  std::set < int > listP;
  for (int c = 0; c < centers_.rows; ++c)
  {
    float min       = FLT_MAX;
    int   pointMin  = -1;
    for (int i = 0; i < points_.rows; ++i){
      if (pointAssignedToCenter_[i] == c){
        if (min > distanceMatrix_(i, c)){
          min       = distanceMatrix_(i, c);
          pointMin  = i;
        }
      }
    }
    if (listP.find(pointMin) != listP.end())
      listP.insert(pointMin);
    if (pointMin != -1){
      for (int j = 0; j < points.cols; ++j){
        centers_(c, j) = points(pointMin, j);
      }
    }
  }
  std::cout << "Cleaning done. Clusters are ready.\n";
}


 //////////////////////////////////////////////////////////////////////////////// 
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
		int indexPoint = ran_.uniform(0, points_.rows);
		while ( listRandom.find(indexPoint) != listRandom.end() )
			indexPoint = ran_.uniform(0, points_.rows);
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
		
		ExecutorService es = Executors.newFixedThreadPool(nbThread_);

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
								
								//Alterado: Verifica se o centro é NaN pois, se não verificar, na tranformação para int o NaN vira Zero e afeta no próximo if
								if(Double.isNaN(center[d])) {
									distance = Double.NaN;
									break; //break pois se o primeiro valor do centro ja é NaN o resto também é NaN
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

////////////////////////////////////////////////////////////////////////////////
//suplmentary function to run the threads 
std::mutex val_mutex;
void KMedians::computeDistanceFor(int ini, int fin)
{
	for (int i = ini; i < fin; ++i)
	{
		for (int c = 0; c < centers_.rows; ++c)
		{
			//if center didn't move, continue
			if(!hasMoved_[c])
				continue;
				
			//if hasn't points in it, continue;
			if(populationInCluster_[c] == 0)
			{
				val_mutex.lock();
				distanceMatrix_(i,c) = FLT_MAX;
				val_mutex.unlock();
				continue;
			}
			
			auto	p		= points_.row(i);
			auto	center	= centers_.row(c);
			double	distance = 0;
			for (int d = 0; d < dimension_; ++d)
			{
				//Alterado: Verifica se o centro é NaN pois, se não verificar, na tranformação para int o NaN vira Zero e afeta no próximo if
				if( center(0,d) == NAN) {
					distance = NAN;
					break; //break pois se o primeiro valor do centro ja é NaN o resto também é NaN
				}
				else 
					distance += kmed_hammingDistanceInt( p(0,d) & 0xff,  (int)center(0,d) & 0xff );	
			}
		}
	}
}
//main computed distace, this function controls the threads 
void KMedians::computeDistance()
{
	int range	= points_.rows / nbThread_,
		n		= 0;
	std::vector<std::thread> vt;
	for (; n < points_.rows; n += range)
		vt.push_back(std::thread ( &KMedians::computeDistanceFor , this, n, n + range) );
	vt.push_back(std::thread ( &KMedians::computeDistanceFor , this, n, points_.rows) );
	for (auto & it: vt){
		it.join();
	}
}
////////////////////////////////////////////////////////////////////////////////
//Este é para a média
/*	private void computeCenters() {
		for(int c = 0 ; c < centers.length; c++)
		{
			if(!hasMoved[c])
				continue;
			if(populationInCluster[c] == 0)
			{
				Arrays.fill(centers[c], Double.NaN);
				continue;
			}
			
			int nbp = 0;
			Arrays.fill(centers[c], 0);
			for(int i = 0 ; i < points.length; i++)
			{
				if(pointAssignedToCenter[i] == c)
				{
					for(int d = 0 ; d < dimension; d++)
						centers[c][d] += points[i][d];
					nbp++;
				}
			}

			if(nbp > 0)
				for(int d = 0 ; d < dimension; d++)
					centers[c][d] = (int)(centers[c][d] / nbp); //Alterado para int truncando para zero assim como no Bag OF Binary Descriptors
			else
				Arrays.fill(centers[c], Double.NaN);
		}
	}
	*/
void KMedians::computeCenters()
{
	for (int c = 0; c < centers_.rows; ++c)
	{
		if (!hasMoved_[c])	continue;
		if (populationInCluster_[c] == 0){
			fillMatRow<float>(centers_, c, NAN);
			continue;
		}
		int nbp = 0;
		fillMatRow<float>(centers_, c, 0);
		for (int i = 0; i < points_.rows; ++i){
			if (pointAssignedToCenter_[i] == c){
				for (int d = 0; d < dimension_; ++d)
					centers_(c,d) += points_(i,d);
				++nbp;
			}
		}
		if (nbp > 0)
		{
			for (int d = 0; d < dimension_; ++d)
				centers_(c, d) = (int)(centers_(c, d) / nbp); //Alterado para int truncando para zero assim como no Bag OF Binary Descriptors
		}
		else
			fillMatRow<float>(centers_, c, NAN);
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
/*
private void computeCentersMedian() {
		for (int c = 0; c < centers.length; c++) {

			if (!hasMoved[c])
				continue;
			if (populationInCluster[c] == 0) {
				Arrays.fill(centers[c], Double.NaN);
				continue;
			}
			
			ArrayList<PairIndexDescriptor<Integer, int[]>> listPointsInClusterC = new ArrayList<PairIndexDescriptor<Integer, int[]>>();
			
			//seleciona todos os pontos do clusters C
			for (int i = 0; i < points.length; i++)
				if (pointAssignedToCenter[i] == c){
					PairIndexDescriptor<Integer, int[]> indexDesc = new PairIndexDescriptor<Integer, int[]>(i, points[i]);
					listPointsInClusterC.add(indexDesc);
				}
			
			//calcula a distância entre todos os pontos do cluster C
			//int distanceMatrixCluster[][] = new int[listPointsInClusterC.size()][listPointsInClusterC.size()];
			int distanceSums[] = new int[listPointsInClusterC.size()];
	
			for(int i = 0; i < listPointsInClusterC.size(); i++)
			{
				//distanceMatrixCluster[i][i] = 0;
				for(int j = i+1; j < listPointsInClusterC.size(); j++)
				{
					int distance = 0;
					for(int d = 0; d < dimension; d++)
					{
						int descI[] = listPointsInClusterC.get(i).Descriptor;
						int descJ[] = listPointsInClusterC.get(j).Descriptor;
						distance += hammingDistance( (((int)descI[d]) & 0xff), (((int)descJ[d]) & 0xff) );
					}
					
					//distanceMatrixCluster[i][j] = distance;
					//distanceMatrixCluster[j][i] = distance;
					
					//distanciâs acumuladas
					distanceSums[i] += distance;
					distanceSums[j] += distance;
				}
			}
			
			int newCenter = -1;
			int shortDist = Integer.MAX_VALUE;
			//Seleciona a mediana do cluster C:
			//É o ponto que minimiza a soma das distâncias aos demais elementos do mesmo cluster
			for(int i = 0; i < distanceSums.length; i++)
			{
				if(distanceSums[i] < shortDist)
				{
					shortDist = distanceSums[i];
					newCenter = i;
				}
			}
			
			int newCenterDesc[] = listPointsInClusterC.get(newCenter).Descriptor;
			for (int j = 0; j < newCenterDesc.length; j++) {
				centers[c][j] = newCenterDesc[j];
			}
			
		}
	}
*/

void KMedians::computeCentersMedian()
{
	for (int c = 0; c < centers_.rows; ++c){
		if (!hasMoved_[c])	continue;
		if (populationInCluster_[c] == 0){
			fillMatRow<float>(centers_, c, NAN);
			continue;
		}
		//ArrayList<PairIndexDescriptor<Integer, int[]>> listPointsInClusterC = new ArrayList<PairIndexDescriptor<Integer, int[]>>();

	}
}

 ////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////Complementary functions/////////////////////////////////////
//simple binary hamming distance between two integers
int kmed_hammingDistanceInt(int x, int y){
	int dist = 0;
	for ( int val = (x ^ y); val; ++dist ) val &= val - 1;		
	return dist;
}
//fill matrix row
template<class t>
void fillMatRow(cv::Mat_<t> & src, int row, t val)
{
	if (row < src.rows)
		for (int i = 0; i < src.cols; ++i)
			src(row, i) = val;
}
