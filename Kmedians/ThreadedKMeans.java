package clustering;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.Collections;

public class ThreadedKMeans implements Serializable, CodeBook
{
	
	private static final long serialVersionUID = -253576488350574703L;

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

	public ThreadedKMeans(int[][] points, int nbCenters, int maxIterations, int nbThreads, Random r)
	{
		this.points = points;
		this.nbCenters = nbCenters;
		this.maxIterations = maxIterations;
		this.nbThread = nbThreads;
		this.ran = r;
		
		//1. init centers
		initCenters();
		
		//2. starting kmeans
		System.out.println("Starting clustering with "+nbThreads+" threads.");
		int nbMoves = 0;
		int iteration = 1;
		do
		{				
			computeCentersMedian();
			
			//compute distance
			computeDistance();
			
			//make assignement
			nbMoves = makeAssignment();
			
			//number of moving centers
			int nbmc = 0;
			for(int c = 0 ; c < centers.length; c++)
				if(hasMoved[c])
					nbmc++;
			
			System.out.println("iteration "+iteration+", "+nbMoves+" points moved, "+nbmc+" centers moved.");
			iteration++;
		}
		while(nbMoves != 0 && iteration < maxIterations);
		System.out.println("Clustering done, cleaning...");
		
		//store variences also
		Arrays.fill(meanDistance, 0);
		for(int c = 0 ; c < centers.length; c++)
		{
			int nbp = 0;
			for(int i = 0; i < points.length; i++)
			{
				if(pointAssignedToCenter[i] == c)
				{
					meanDistance[c] += distanceMatrix[i][c];
					nbp ++;
				}
			}
			if(nbp > 0)
				meanDistance[c] /= (double) nbp;
		}
		
		//cleaning empty clusters
		ArrayList<double[]> listOfCenters = new ArrayList<double[]>();
		ArrayList<Double> listOfMeanDist = new ArrayList<Double>();
		ArrayList<Integer> listOfPopulation = new ArrayList<Integer>();
		for(int c = 0 ; c < centers.length; c++)
		{
			if(populationInCluster[c] > 0)
			{
				listOfCenters.add(centers[c]);
				listOfMeanDist.add(meanDistance[c]);
				listOfPopulation.add(populationInCluster[c]);
			}
		}
		
		centers = new double[listOfCenters.size()][];
		meanDistance = new double[listOfCenters.size()];
		populationInCluster = new int[listOfCenters.size()];
		for(int c = 0 ; c < centers.length; c++)
		{
			centers[c] = listOfCenters.get(c);
			meanDistance[c] = listOfMeanDist.get(c);
			populationInCluster[c] = listOfPopulation.get(c);
		}
		
		ArrayList<Integer> listP = new ArrayList<Integer>();
		for(int c = 0 ; c < centers.length; c++)
		{
			double min = Double.POSITIVE_INFINITY;
			int pointMin = -1;
			for(int i = 0; i < points.length; i++) {
				if(pointAssignedToCenter[i] == c) {
					if (min > distanceMatrix[i][c]) {
						min = distanceMatrix[i][c];
						pointMin = i;
					}
				}
			}
			if (!listP.contains(pointMin))
				listP.add(pointMin);
			if (pointMin!=-1) {
				for (int j = 0 ; j < points[0].length; j++)
					centers[c][j] = points[pointMin][j];
			}
		}
		
		System.out.println("Cleaning done. Clusters are ready.");
	}
	
	//Para computar com mediana
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
		
	//Este é para a média
	private void computeCenters() {
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

	private void computeDistance() {
		
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
	
	private int makeAssignment() {

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
				
				/*
				ArrayList<Integer> a = pointsInCluster.get(indexMin);
				a.add(i);
				pointsInCluster.set(indexMin, a);
				*/
				
			}
		}	
		return nbm;
	}

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
			
			/*
			ArrayList<Integer> p =  new ArrayList<Integer>();
			p.add(indexPoint);
			pointsInCluster.set(i, p);
			*/
			
			populationInCluster[i]++;
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
	
	/*private boolean isCenter(int x)
	{
		for(int i = 0; i < center_Point.length; i++)
		{
			if(center_Point[i] == x)
				return true;
		}
		return false;
	}*/
	
	private int hammingDistance(int x, int y) 
	{
		int dist = 0;
		int val = x ^ y; // XOR

		// Count the number of set bits
		while (val != 0) {
			++dist;
			val &= val - 1;
		}

		return dist;
	}

	/**
	 * gets the centers of the clusters
	 * @return
	 */
	public double[][] getClusters()
	{
		return centers;
	}
	
	/**
	 * get the mean distance of each each cluster regarding its points
	 * @return
	 */
	public double[][] getCovarianceMatrices() {
//		return meanDistance;
		return null;
	}

	//@Override
	public double[] getClusterRadii() {
		return meanDistance;
	}
	
	//@Override
	public int getNumberOfCenters() {
		return nbCenters;
	}
	
	/**
	 * get the number of points in each cluster
	 * @return
	 */
	public int[] getPopulationInCluster() {
		return populationInCluster;
	}
}