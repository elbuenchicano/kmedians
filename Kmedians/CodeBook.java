package clustering;

/**
 * @author picard
 * 
 * Interface for all clustering utilities
 */
public interface CodeBook {
	
	/**
	 * get the centers of all clusters
	 * @return
	 */
	public double[][] getClusters();
	
	/**
	 * get covariance matrices of all clusters
	 * @return
	 */
	public double[][] getCovarianceMatrices();
	
	/**
	 * get radius (mean of distances) of all clusters
	 * @return
	 */
	public double[] getClusterRadii();
	
	/**
	 * get the number of items falling in each cluster
	 * @return
	 */
	public int[] getPopulationInCluster();
	
	/**
	 * get the number of centers (visual words)
	 * @return
	 */
	public int getNumberOfCenters();

}
