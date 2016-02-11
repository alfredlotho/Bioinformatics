import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;

public class NeighborJoining {
	/**
	 * 
	 * @param distMatrix - an additive distance matrix
	 * @return a distance matrix (D*) where every element is derived using the formula:
	 * 	D*[i][j] = (n*2)*D[i][j] - TotalDistance(i) - TotalDistance(j)
	 */
	public static double[][] ConstructNeighborJoiningMatrix(double[][] distMatrix, int actualSize) {
		int size = distMatrix[0].length;
		double[][] dStar = new double[size][size];
		
		double[] totalDistances = new double[size];
		
		for (int i = 0; i < size; i++) {
			totalDistances[i] = TotalDistance(distMatrix[i]);
		}
		
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (distMatrix[i][j] < 0) {
					dStar[i][j] = -1; //set to empty for sparse matrix (a matrix whose rows or columns have been marked as empty by setting to -1)
				} else if (distMatrix[i][j] == 0) {
					dStar[i][j] = 0;
				} else {
					dStar[i][j] = ((actualSize-2) * distMatrix[i][j]) - totalDistances[i] - totalDistances[j];
				}
			}
		}
		
		return dStar;
	}
	
	/**
	 * 
	 * @param row - a row in the distance matrix
	 * @return the sum of all elements in the row
	 */
	private static double TotalDistance(double[] row) {
		double totalDist = 0;
		for (int i = 0; i < row.length; i++) {
			if (row[i] > 0) { //0 is the identity line, -1 is an empty cell in a sparse matrix
				totalDist += row[i];
			}
		}
		return totalDist;
	}
	
	/**
	 * 
	 * @param dStar - the neighbor joining distance matrix
	 * @return the row and column index of the element with the smallest value
	 */
	private static Cluster FindMinNonDiagonalElement(double[][] dStar) {
		double minVal = Double.MAX_VALUE;
		int minI = 0;
		int minJ = 0;
		int size = dStar[0].length;
		
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < size; j++) {
				if (dStar[i][j] < minVal && i != j 
						&& dStar[i][j] != -1) {
					minVal = dStar[i][j];
					minI = i;
					minJ = j;
				}
			}
		}
		
		return new Cluster(minI, minJ, minVal);
	}
	
	/**
	 * 
	 * @param dStar - a 2x2 sparse neighbor joining matrix containing -1 for empty elements
	 * @return the index of the remaining postive element in the sparse matrix
	 */
	private static Cluster FindRemainingPositiveElement(double[][] dStar) {
		int size = dStar[0].length;
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < i; j++) {
				if (dStar[i][j] > 0) {
					return new Cluster(i, j, dStar[i][j]);
				}
			}
		}
		return null;
	}
	
	static NeighborTree neighborTree;
	public static NeighborTree ExecuteJoiningAlgorithm(double[][] d, int n, int m) {
		if (n == 2) {
			Cluster c = FindRemainingPositiveElement(d);
			return new NeighborTree(c.nodeA, c.nodeB, c.weight);
		}
		double[][] dStar = ConstructNeighborJoiningMatrix(d, n);
		Cluster c = FindMinNonDiagonalElement(dStar);
		double delta = (TotalDistance(d[c.nodeA]) - TotalDistance(d[c.nodeB])) / (n-2);
		double limbLength1 = (d[c.nodeA][c.nodeB] + delta)/2;
		double limbLength2 = (d[c.nodeA][c.nodeB] - delta)/2;
		double[][] dReduced = CombineCluster(d, c.nodeA, c.nodeB);
		neighborTree = ExecuteJoiningAlgorithm(dReduced, n-1, m+1);
		neighborTree.AddNode(m);
		neighborTree.AddEdge(m, c.nodeA, limbLength1);
		neighborTree.AddEdge(m, c.nodeB, limbLength2);
		return neighborTree;
	}
	
	public static double[][] CombineCluster(double[][] origMatrix, int i, int j) {
		int origSize = origMatrix[0].length;
		double[][] expandedMatrix = new double[origSize+1][origSize+1];
				
		for (int k = 0; k < origSize; k++) {
			if (k == i || k == j) {
				expandedMatrix[origSize][k] = -1;
				expandedMatrix[k][origSize] = -1;
			} else {
				expandedMatrix[origSize][k] = (origMatrix[i][k] + origMatrix[j][k] - origMatrix[i][j]) / 2;
				expandedMatrix[k][origSize] = expandedMatrix[origSize][k];
			}
		}
		
		for (int r = 0; r < origSize; r++) {
			for (int c = 0; c < origSize; c++) {
				if (r == i || c == i || r == j || c == j) {
					expandedMatrix[r][c] = -1;
				} else {
					expandedMatrix[r][c] = origMatrix[r][c];
				}
			}
		}
		expandedMatrix[origSize][origSize] = 0;
		
		return expandedMatrix;
	}
}

class Cluster {
	int nodeA;
	int nodeB;
	double weight;
	
	public Cluster(int nodeA, int nodeB, double weight) {
		this.nodeA = nodeA;
		this.nodeB = nodeB;
		this.weight = weight;
	}
}


class NeighborNode {
	int index;
	
	public NeighborNode(int index) {
		this.index = index;
	}
}

class NeighborEdge {
	int nodeA;
	int nodeB;
	double weight;
	
	public NeighborEdge(int nodeA, int nodeB, double weight) {
		this.nodeA = nodeA;
		this.nodeB = nodeB;
		this.weight = weight;
	}
	
	@Override
	public String toString() {
		DecimalFormat df = new DecimalFormat("#.000");
		return nodeA + BioinformaticsCommon.NODE_SEPARATOR +nodeB +":" +df.format(weight);
	}
}

class NeighborTree {
	List<NeighborNode> nodeList;
	List<NeighborEdge> edgeList;
	
	public NeighborTree(int nodeA, int nodeB, double weight) {
		nodeList = new ArrayList<NeighborNode>();
		edgeList = new ArrayList<NeighborEdge>();
		AddNode(nodeA);
		AddNode(nodeB);
		AddEdge(nodeA, nodeB, weight);
	}
	
	public void AddNode(int nodeIndex) {
		if (!nodeList.contains(nodeIndex))
			nodeList.add(new NeighborNode(nodeIndex));
	}
	
	public void AddEdge(int nodeA, int nodeB, double weight) {
		edgeList.add(new NeighborEdge(nodeA, nodeB, weight));
		edgeList.add(new NeighborEdge(nodeB, nodeA, weight));
	}
	
	public void PrintTree() {
		Comparator<NeighborEdge> byNodeA = new Comparator<NeighborEdge>() {
	        @Override
	        public int compare(NeighborEdge x, NeighborEdge y) {
	        	return ((Integer)x.nodeA).compareTo(y.nodeA);
	        }
	    };
		Collections.sort(edgeList, byNodeA);
		for (int i = 0; i < edgeList.size(); i++) {
			System.out.println(edgeList.get(i));
		}
	}
}