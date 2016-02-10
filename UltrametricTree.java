import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

public class UltrametricTree {
	List<UltrametricNode> nodeList;
	List<UltrametricEdge> edgeList;
	
	public UltrametricTree() {
		nodeList = new ArrayList<UltrametricNode>();
		edgeList = new ArrayList<UltrametricEdge>();
	}
	
	public void UPGMA(int leafCount, double[][] distMatrix) {
		for (int i = 0; i < leafCount; i++) {
			nodeList.add(new UltrametricNode(0));
		}
		double[][] reducedMatrix = distMatrix;
		for (int i = leafCount; i <= 2*(leafCount-1); i++) {
			reducedMatrix = FindClosestCluster(reducedMatrix);
		}		
	}
	
	public void PrintEdges() {
		Comparator<UltrametricEdge> byNodeA = new Comparator<UltrametricEdge>() {
	        @Override
	        public int compare(UltrametricEdge x, UltrametricEdge y) {
	        	return ((Integer)x.nodeA).compareTo(y.nodeA);
	        }
	    };
		Collections.sort(edgeList, byNodeA);
		for (int i = 0; i < edgeList.size(); i++) {
			System.out.println(edgeList.get(i));
		}
	}
	
	public double[][] FindClosestCluster(double[][] distMatrix) {
		int minI = -1;
		int minJ = -1;
		double minDist = Double.MAX_VALUE;
		for (int i = 0; i < distMatrix[0].length; i++) {
			for (int j = 0; j < i; j++) {
				if (distMatrix[i][j] < minDist && distMatrix[i][j] != -1) {
					minDist = distMatrix[i][j];
					minI = i;
					minJ = j;
				}
			}
		}
		nodeList.add(new UltrametricNode(minDist/2.0));
		int parentNode = nodeList.size()-1; //the last added node
		nodeList.get(parentNode).effectiveNodeCount = nodeList.get(minI).effectiveNodeCount + nodeList.get(minJ).effectiveNodeCount;
		double weight1 = nodeList.get(parentNode).age - nodeList.get(minI).age;
		edgeList.add(new UltrametricEdge(parentNode, minI, weight1));
		edgeList.add(new UltrametricEdge(minI, parentNode, weight1));
		double weight2 = nodeList.get(parentNode).age - nodeList.get(minJ).age;
		edgeList.add(new UltrametricEdge(parentNode, minJ, weight2));
		edgeList.add(new UltrametricEdge(minJ, parentNode, weight2));
		return CombineCluster(distMatrix, minI, minJ);
	}
	
	public double[][] CombineCluster(double[][] origMatrix, int minI, int minJ) {
		int origSize = origMatrix[0].length;
		double[][] expandedMatrix = new double[origSize+1][origSize+1];
				
		for (int c = 0; c < origSize; c++) {
			if (c == minI || c == minJ) {
				expandedMatrix[origSize][c] = -1;
				expandedMatrix[c][origSize] = -1;
			} else {
				int iCount = nodeList.get(minI).effectiveNodeCount;
				int jCount = nodeList.get(minJ).effectiveNodeCount;
				expandedMatrix[origSize][c] = (origMatrix[minI][c]*iCount + origMatrix[minJ][c]*jCount)/(iCount + jCount);
				expandedMatrix[c][origSize] = expandedMatrix[origSize][c];
			}
		}
		
		for (int r = 0; r < origSize; r++) {
			for (int c = 0; c < origSize; c++) {
				if (r == minI || c == minI || r == minJ || c == minJ) {
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

class UltrametricNode {
	double age;
	int effectiveNodeCount = 1; //increase node count if this node is a result of merging two nodes
	
	public UltrametricNode(double age) {
		this.age = age;
	}
}

class UltrametricEdge {
	int nodeA;
	int nodeB;
	double weight;
	
	public UltrametricEdge(int nodeA, int nodeB, double weight) {
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
