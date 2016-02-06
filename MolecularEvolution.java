import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

public class MolecularEvolution {

	public static void main(String[] args) {
		// common
		String path = "C:\\Users\\Timothy\\Desktop\\temp\\bioinformatics\\dataset_a_b.txt";
		List<String> lines = BioinformaticsCommon.readFile(path);
		
		// test for AdditivePhylogeny
		int n = Integer.parseInt(lines.get(0));
		origN = n;
		int[][] dist = ConstructDistanceMatrixFromFile(lines, 1);
		newNodeLabel = n*2-3;
		PhylogenicTree tree = AdditivePhylogeny(dist, n-1);
		List<Map.Entry<String, Edge>> list = tree.GetSortedEdges();
		for (int s = 0; s < list.size(); s++) {
			Map.Entry me2 = (Map.Entry)list.get(s);
			System.out.print(me2.getKey() + ":" +((Edge)me2.getValue()).weight +"\r\n");
		}
	}
	
	/**
	 * generates an adjacency matrix from a list, each row in the format a->b:c where a and b
	 * are the connected nodes and c represents the weight
	 */
	private static int[][] ConvertAdjacencyListToMatrix(List<String> lines, int leafCount) {
		int maxNodeCount = lines.size() - 1 - leafCount;
		int[][] adjacencyMatrix = new int[maxNodeCount][maxNodeCount];
		
		for (int lineIndex = 1; lineIndex < lines.size(); lineIndex++) {
			StringTokenizer st = new StringTokenizer(lines.get(lineIndex), ":");
			StringTokenizer st2 = new StringTokenizer(st.nextToken(), PhylogenicTree.NODE_SEPARATOR);
			int weight = Integer.parseInt(st.nextToken());
			int leftNode = Integer.parseInt(st2.nextToken());
			int rightNode = Integer.parseInt(st2.nextToken());
			adjacencyMatrix[leftNode][rightNode] = weight;
		}
		return adjacencyMatrix;
	}
	
	/**
	 * Uses the floyd-warshall algorithm to compute the shortest distance between nodes given a graph
	 * Distances Between Leaves Problem: Compute the distances between leaves in a weighted tree.
     * Input:  An integer n followed by the adjacency list of a weighted tree with n leaves.
     * Output: An n x n matrix (di,j), where di,j is the length of the path between leaves i and j.
     * Sample Input:
	 *	4
	 *	0->4:11
	 *	1->4:2
	 *	2->5:6
	 *	3->5:7
	 *	4->0:11
	 *	4->1:2
	 *	4->5:4
	 *	5->4:4
	 *	5->3:7
	 *	5->2:6
	 *	Sample Output:
	 *	0	13	21	22
	 *	13	0	12	13
	 *	21	12	0	13
	 *	22	13	13	0
	 */
	private static void DistanceBetweenLeaves(List<String> lines) {
		int leafCount = Integer.parseInt(lines.get(0));
		int[][] adjacencyMatrix = ConvertAdjacencyListToMatrix(lines, leafCount);
		int size = adjacencyMatrix[0].length;
		int MAX = 9999; //do not use Integer.MAX because it fails when it is compared
		
		int[][] dist = new int[size][size];
		for (int k = 0; k < size; k++) {
			for (int i = 0; i < size; i++) {
				for (int j = 0; j < size; j++) {
					if (i == j) {
						dist[i][j] = 0;
					} else if (dist[i][j] == 0) {
						if (adjacencyMatrix[i][j] == 0) {
							dist[i][j] = MAX;
						} else {
							dist[i][j] = adjacencyMatrix[i][j];
						}
					} else if (dist[i][j] > dist[i][k] + dist[k][j]) {
						dist[i][j] = dist[i][k] + dist[k][j];
					}
				}
			}
		}
		
		for (int i = 0; i < leafCount; i++) {
			for (int j = 0; j < leafCount; j++) {
				System.out.print(dist[i][j] + " ");
			}
			System.out.println();
		}
	}
	
	/**
	 * Used for converting the file input to a 2d distance matrix
	 * @param lines - line by line parsed string in the following format:
	 * 	line 0: the number of leaves in the graph; also corresponds to N, the size of the distance matrix
	 * 	line 1: the index of the leaf we are getting a limb length for
	 *  line 2-n: a space separated N x N matrix containing the distances between each of the nodes in a graph
	 * @param startWithLineIndex - the 0 based index of the line where the matrix starts
	 * @return
	 */
	private static int[][] ConstructDistanceMatrixFromFile(List<String> lines, int startWithLineIndex) {
		int leafCount = Integer.parseInt(lines.get(0));
		int[][] dist = new int[leafCount][leafCount];
		
		for (int lineIndex = startWithLineIndex; lineIndex < lines.size(); lineIndex++) {
			StringTokenizer st = new StringTokenizer(lines.get(lineIndex));
			for (int n = 0; n < leafCount; n++) {
				dist[lineIndex-startWithLineIndex][n] = Integer.parseInt(st.nextToken());
			}
		}
		
		return dist;
	}
	
	/**
	 * Used for computing the length of a limb (the distance from a leaf j to its parent)
	 * @param dist - the distance matrix
	 * @param j - the random leaf to be used as reference point 
	 * @return the distance between leaf j and its parent node
	 */
	private static int LimbLength(int[][] dist, int j) {
		int limbLength = Integer.MAX_VALUE;
		int size = dist[0].length;
		
		for (int i = 0; i < size; i++) {
			for (int k = 0; k < size; k++) {
				if (i != j && j != k) {
					limbLength = Math.min(limbLength, (dist[i][j] + dist[j][k] - dist[i][k]) / 2);
				}
			}
		}
		
		return limbLength;
	}
	
	

	/**
	 * 
	 * @param dist - the distance matrix
	 * @param n - equivalent to leafCount - 1 (this is the last node)
	 * @return
	 */
	static int[][] dTrimmed;
	static int[][] dBald;
	static int newNodeLabel; //initialize to leafCount before calling function
	static PhylogenicTree tree;
	static int origN;
	private static PhylogenicTree AdditivePhylogeny(int[][] dist, int n) {
		if (n == 1) {
			return new PhylogenicTree(dist[0][1], origN);
		}
		int limbLength = LimbLength(dist, n);
		for (int j = 0; j < n; j++) {
			dist[j][n] -= limbLength;
			dist[n][j] = dist[j][n];
		}
		
		// find three leaves (i, n, k) 
		int left = 0;
		int right = 0;
		int x = 0;
		int y = 0;
		for (int i = 0; i < n; i++) {
			for (int k = 0; k < i; k++) {
				if (dist[i][k] == dist[i][n] + dist[n][k]) {
					x = dist[i][n];
					y = dist[i][k];
					left = i;
					right = k;
					break;
				}
			}
		}
		dTrimmed = ReduceMatrix(dist, 1);
		tree = AdditivePhylogeny(dTrimmed, n-1);
		tree.BreakEdge(left, right, newNodeLabel, x, y, n, limbLength);
		newNodeLabel--;
		return tree;
	}

	
	/**
	 * 
	 * @param dist - a distance matrix of size N by N
	 * @return the matrix of size N-1 by N-1 (removing the last row and last column)
	 */
	private static int[][] ReduceMatrix(int[][] dist, int count) {
		int reducedLength = dist[0].length - count;
		int[][] reducedMatrix = new int[reducedLength][reducedLength];
		for (int i = 0; i < reducedLength; i++) {
			System.arraycopy(dist[i], 0, reducedMatrix[i], 0, dist[i].length-count);
		}
		return reducedMatrix;
	}
}