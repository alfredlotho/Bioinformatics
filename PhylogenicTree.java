import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

public class PhylogenicTree {
	Map<Integer, Node> nodeList = new HashMap<Integer, Node>();
	Map<String, Edge> edgeList = new HashMap<String, Edge>();
	int[][] pathMatrix;
	int origN = 0;
	
	public PhylogenicTree(int initialWeight, int origN) {
		this.origN = origN;
		int emSize = (origN * 2) - 2;
		pathMatrix = new int[emSize][emSize];
		AddToPathMatrix(0, 1, initialWeight);
		nodeList.put(0, new Node(0, 1));
		nodeList.put(1, new Node(1, 0));
		edgeList.put("0" +BioinformaticsCommon.NODE_SEPARATOR +"1", new Edge(nodeList.get(0), nodeList.get(1), initialWeight));
		edgeList.put("1" +BioinformaticsCommon.NODE_SEPARATOR +"0", new Edge(nodeList.get(1), nodeList.get(0), initialWeight));
	}
	
	public void AddToPathMatrix(int a, int b, int weight) {
		pathMatrix[a][b] = weight;
		pathMatrix[b][a] = weight;
	}
	
	public String FindMiddleNode(int a, int b, int weight, int rightBreakWeight) {
		for (int i = 0; i < pathMatrix[0].length; i++) {
			if (pathMatrix[a][i] < weight && pathMatrix[i][b] != 0 && i >= origN) {
				for (int j = 0; j < pathMatrix[0].length; j++) {
					if (pathMatrix[i][j] != 0 && (pathMatrix[j][b] != 0 || j == b) &&
						pathMatrix[i][j] + pathMatrix[j][b] == pathMatrix[i][b] &&
						pathMatrix[a][j] + pathMatrix[j][b] == pathMatrix[a][b] &&
						pathMatrix[a][j] > weight) {		
							if (!edgeList.containsKey(i +BioinformaticsCommon.NODE_SEPARATOR +j)) {
								continue;
							}
							int nodeWeight = weight - pathMatrix[a][i];
							return i +" " +j +" " +nodeWeight;
						}
				}
			}
		}
		return "";
	}
	
	public List GetSortedEdges() {
		Comparator<Map.Entry<String, Edge>> byNodeA = new Comparator<Map.Entry<String, Edge>>() {
	        @Override
	        public int compare(Map.Entry<String, Edge> left, Map.Entry<String, Edge> right) {
	        	Integer leftComparator = Integer.parseInt(new StringTokenizer(left.getKey(), BioinformaticsCommon.NODE_SEPARATOR).nextToken());
	        	Integer rightComparator = Integer.parseInt(new StringTokenizer(right.getKey(), BioinformaticsCommon.NODE_SEPARATOR).nextToken());
	            return leftComparator.compareTo(rightComparator);
	        }
	    };
	    
		List<Map.Entry<String, Edge>> list = new ArrayList<Map.Entry<String, Edge>>(edgeList.entrySet());
		Collections.sort(list, byNodeA);
		return list;
	}
	
	public void BreakEdge(int leftNodeIndex, int rightNodeIndex, int newNodeId, int leftBreakWeight, int totalWeight, int baldIndex, int limbLength) {
		int nodeAIndex = leftNodeIndex;
		int nodeBIndex = nodeList.get(nodeAIndex).parentNode;
		String edgeId = nodeAIndex +BioinformaticsCommon.NODE_SEPARATOR +nodeBIndex;
		if (edgeList.get(edgeId).weight < leftBreakWeight) {
			StringTokenizer st = new StringTokenizer(FindMiddleNode(leftNodeIndex, rightNodeIndex, leftBreakWeight, totalWeight-leftBreakWeight));
			nodeAIndex = Integer.parseInt(st.nextToken());
			nodeBIndex = Integer.parseInt(st.nextToken());
			leftBreakWeight = Integer.parseInt(st.nextToken());
			edgeId = nodeAIndex +BioinformaticsCommon.NODE_SEPARATOR +nodeBIndex;
		}
		Node newParentNode = new Node(newNodeId, nodeBIndex);
		// break the edge by dividing the distance between the node to its parent by value of leftBreakWeight
		nodeList.put(newNodeId, newParentNode);
		edgeList.put(nodeAIndex+BioinformaticsCommon.NODE_SEPARATOR+newNodeId, new Edge(nodeList.get(nodeAIndex), newParentNode, leftBreakWeight));
		edgeList.put(newNodeId+BioinformaticsCommon.NODE_SEPARATOR+nodeAIndex, new Edge(newParentNode, nodeList.get(nodeAIndex), leftBreakWeight));
		nodeList.get(nodeAIndex).parentNode = newNodeId;
		AddToPathMatrix(nodeAIndex, newNodeId, leftBreakWeight);
		// connect new internal node to the right part of the broken edge
		int rightBreakWeight = edgeList.get(edgeId).weight - leftBreakWeight;
		edgeList.put(nodeBIndex+BioinformaticsCommon.NODE_SEPARATOR+newNodeId, new Edge(nodeList.get(nodeBIndex), newParentNode, rightBreakWeight));
		edgeList.put(newNodeId+BioinformaticsCommon.NODE_SEPARATOR+nodeBIndex, new Edge(newParentNode, nodeList.get(nodeBIndex), rightBreakWeight));
		nodeList.get(nodeBIndex).parentNode = newNodeId;
		AddToPathMatrix(nodeBIndex, newNodeId, rightBreakWeight);
		// connect new internal node to bald limb node
		nodeList.put(baldIndex, new Node(baldIndex, newNodeId));
		edgeList.put(baldIndex+BioinformaticsCommon.NODE_SEPARATOR+newNodeId, new Edge(nodeList.get(baldIndex), newParentNode, limbLength));
		edgeList.put(newNodeId+BioinformaticsCommon.NODE_SEPARATOR+baldIndex, new Edge(newParentNode, nodeList.get(baldIndex), limbLength));
		AddToPathMatrix(baldIndex, newNodeId, limbLength);
		// remove original edge that was broken into 3 parts
		edgeList.remove(edgeId);
		edgeList.remove(nodeBIndex+BioinformaticsCommon.NODE_SEPARATOR+nodeAIndex);
		
		// recalculate path matrix
		for (int z : nodeList.keySet()) {
			if (z != baldIndex && z != newNodeId) {
				if (pathMatrix[nodeAIndex][z] != 0 && pathMatrix[nodeAIndex][z] < pathMatrix[nodeBIndex][z]) {
						pathMatrix[newNodeId][z] = pathMatrix[newNodeId][nodeAIndex] + pathMatrix[nodeAIndex][z];
						pathMatrix[z][newNodeId] = pathMatrix[newNodeId][z];
				}
				if (pathMatrix[nodeBIndex][z] != 0 && pathMatrix[nodeBIndex][z] < pathMatrix[nodeAIndex][z]) {
					pathMatrix[newNodeId][z] = pathMatrix[newNodeId][nodeBIndex] + pathMatrix[nodeBIndex][z];
					pathMatrix[z][newNodeId] = pathMatrix[newNodeId][z];
				}
			}
		}
		for (int z : nodeList.keySet()) {
			if (z != baldIndex && z != newNodeId) {
				pathMatrix[baldIndex][z] = pathMatrix[baldIndex][newNodeId] + pathMatrix[newNodeId][z];
				pathMatrix[z][baldIndex] = pathMatrix[baldIndex][z];
			}
		}
		
		
	}
}

class Node {
	int index;
	int parentNode;
	
	public Node(int index, int parentNode) {
		this.index = index;
		this.parentNode = parentNode;
	}
}

class Edge {
	Node nodeA;
	Node nodeB;
	int weight;
	
	public Edge(Node nodeA, Node nodeB, int weight) {
		this.nodeA = nodeA;
		this.nodeB = nodeB;
		this.weight = weight;
	}
	
}