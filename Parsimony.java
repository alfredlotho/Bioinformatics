import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

public class Parsimony {
	int dnaLength;
	Map<Integer, ParsimonyTree> treeList;
	static final char[] ALPHABET = {'A','C','G','T'};
	static final char NULL_NUCLEOTIDE = '-';
	static final int INFINITY = 99999; //do not use Integer.MAX_VALUE since comparison error occurs
	static boolean isScoreShown = true;
	StringBuffer sb = new StringBuffer();
	
	/**
	 * performs the small parsimony algorithm on each of the character indices of a rooted tree
	 * (i.e. if there are 8 characters on the dna string, this algorithm is performed 8 times and then the
	 * all the results are combined together to form the dna string for each internal node (and root) that 
	 * produces the minimum parsimony score for the input adjacency list)
	 * @param lines - a number on the first line representing the number of leaves in the tree followed by
	 * 					an adjacency list representing the predefined rooted tree structure
	 */
	public void RootedSmallParsimony(List<String> lines) {
		treeList = new HashMap<Integer, ParsimonyTree>();
		int leafCount = Integer.parseInt(lines.get(0));
		ParsimonyTree mainTree = new ParsimonyTree(leafCount);
		for (int i = 1; i < lines.size(); i++) {
			StringTokenizer st = new StringTokenizer(lines.get(i), BioinformaticsCommon.NODE_SEPARATOR);
			int parentIndex = Integer.parseInt(st.nextToken());
			if (i <= leafCount) {
				String childLabel = st.nextToken();
				
				if (dnaLength == 0) {
					dnaLength = childLabel.length();
					treeList.put(0, mainTree);
					for (int j = 1; j <= dnaLength; j++) {
						treeList.put(j, new ParsimonyTree(leafCount));
					}
					
				}
				mainTree.AddNode(i-1, -1, childLabel);
				mainTree.AddNode(parentIndex, i-1, "");
				
				for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
					treeList.get(charIndex+1).AddNode(i-1, -1, childLabel.substring(charIndex, charIndex+1));
					treeList.get(charIndex+1).AddNode(parentIndex, i-1, "");
				}
			} else {
				int childIndex = Integer.parseInt(st.nextToken());
				mainTree.AddNode(parentIndex, childIndex, "");
				for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
					treeList.get(charIndex+1).AddNode(parentIndex, childIndex, "");
				}
			}
		}
		
		for (int k = 1; k <= dnaLength; k++) {
			treeList.get(k).SmallParsimony();
			for (int m = leafCount; m < mainTree.nodeList.size(); m++) {
				mainTree.nodeList.get(m).label += treeList.get(k).nodeList.get(m).label;
			}
		}
		mainTree.PrintEdges(true, true, sb);
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
	}
	
	public ParsimonyTree UnrootedSmallParsimony(List<String> lines) {
		treeList = new HashMap<Integer, ParsimonyTree>();
		int leafCount = Integer.parseInt(lines.get(0));
		ParsimonyTree mainTree = new ParsimonyTree(leafCount);
		int parentIndex;
		
		int leafIndex = -1;
		for (int i = 1; i < lines.size(); i++) {
			StringTokenizer st = new StringTokenizer(lines.get(i), BioinformaticsCommon.NODE_SEPARATOR);
			String leftElem = st.nextToken();
			if (!BioinformaticsCommon.isInteger(leftElem)) {
				continue; //ignore this line if it starts with the dna string since this is a duplicate edge
			}
			parentIndex = Integer.parseInt(leftElem);
			String childLabel = st.nextToken();
			if (!BioinformaticsCommon.isInteger(childLabel)) {
				leafIndex++;
					
				if (dnaLength == 0) {
					dnaLength = childLabel.length();
					treeList.put(0, mainTree);
					for (int j = 1; j <= dnaLength; j++) {
						treeList.put(j, new ParsimonyTree(leafCount));
					}
					
				}
				mainTree.AddNode(leafIndex, -1, childLabel);
				mainTree.AddNode(parentIndex, leafIndex, "");
				
				
				for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
					treeList.get(charIndex+1).AddNode(leafIndex, -1, childLabel.substring(charIndex, charIndex+1));
					treeList.get(charIndex+1).AddNode(parentIndex, leafIndex, "");
				}
			} else { 
				int childIndex = Integer.parseInt(childLabel);
				mainTree.AddNode(parentIndex, childIndex, "");
				for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
					treeList.get(charIndex+1).AddNode(parentIndex, childIndex, "");
				}
			}
		}
						
		// add the final root index that will be removed later
		parentIndex = mainTree.nodeList.size();
		// connect daughter and son node to temporary root node
		InsertRoot(mainTree, parentIndex);
		for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
			InsertRoot(treeList.get(charIndex+1), parentIndex);
		}

		for (int tIndex = 0; tIndex < treeList.size(); tIndex++) {
			ParsimonyTree tree = treeList.get(tIndex);
			for (int pIndex = tree.nodeList.size()-1; pIndex >= leafCount; pIndex--) {
				ParsimonyNode currNode = tree.nodeList.get(pIndex);
				ParsimonyNode a = currNode.connection[0];
				ParsimonyNode b = currNode.connection[1];
				if (a == null)
					; // do nothing
				else if (a.index < leafCount)
					a.connection[2] = currNode;
				else
					tree.AddNode(a.index, currNode.index, "");
				
				if (b == null)
					; // do nothing
				else if (b.index < leafCount)
					b.connection[2] = currNode; 
				else
					tree.AddNode(b.index, currNode.index, "");
				
				if (!currNode.isRoot()) {
					ParsimonyNode c = currNode.connection[2];
					if (c == null)
						; // do nothing
					else if (c.index < leafCount)
						c.connection[2] = currNode; 
					else
						tree.AddNode(c.index, currNode.index, "");
				}
			}
		}
		
		for (int k = 1; k <= dnaLength; k++) {
			treeList.get(k).SmallParsimony();
			for (int m = leafCount; m < mainTree.nodeList.size(); m++) {
				mainTree.nodeList.get(m).label += treeList.get(k).nodeList.get(m).label;
			}
		}
		//mainTree.PrintNodes();
		
		return mainTree;
	}
	
	/**
	 * Insert the root node between a random edge
	 * Edge selection strategy: the last added non-root node and an internal node connected to it
	 * @param tree - the tree where the changes will be reflected
	 * @param parentIndex - the index of the newly created temporary root node
	 */
	private void InsertRoot(ParsimonyTree tree, int parentIndex) {
		int daughter = parentIndex - 1;
		ParsimonyNode daughterNode = tree.nodeList.get(parentIndex - 1);
		int son = -1;
		// find an internal node connected to daughter
		for (int i = 0; i < daughterNode.connection.length; i++) {
			if (!daughterNode.connection[i].isLeaf()) {
				son = daughterNode.connection[i].index;
				break;
			}
		}
		tree.tempRootInsPoint = new TemporaryRootInsertionPoint(daughter, son);
		tree.AddNode(parentIndex, tree.tempRootInsPoint.nodeA, "");
		tree.AddNode(parentIndex, tree.tempRootInsPoint.nodeB, "");
		// assign the temporary root node as one of the connections of both nodes in the root insertion point  
		ReplaceParent(tree.nodeList.get(tree.tempRootInsPoint.nodeA), 
				tree.nodeList.get(tree.tempRootInsPoint.nodeB), 
				tree.nodeList.get(parentIndex),
				true);
		ReplaceParent(tree.nodeList.get(tree.tempRootInsPoint.nodeB), 
				tree.nodeList.get(tree.tempRootInsPoint.nodeA), 
				tree.nodeList.get(parentIndex),
				true);
	}
	
	/**
	 * Input: Two internal nodes a and b specifying an edge e, followed by an adjacency
     * list of an unrooted binary tree.
     * Output: Two adjacency lists representing the nearest neighbors of the tree with
     * respect to e. The adjacency lists are separated with a blank line.
	 */
	public void NearestNeighborsOfTree(List<String> lines) {
		StringTokenizer st = new StringTokenizer(lines.get(0));
		int internalNodeA = Integer.parseInt(st.nextToken());
		int internalNodeB = Integer.parseInt(st.nextToken());
		
		ParsimonyTree tree = new ParsimonyTree(-1);
		for (int i = 1; i < lines.size(); i++) {
			st = new StringTokenizer(lines.get(i), BioinformaticsCommon.NODE_SEPARATOR);
			int left = Integer.parseInt(st.nextToken());
			int right = Integer.parseInt(st.nextToken());
			tree.AddNode(right, -1, "");
			tree.AddNode(left, right, "");
		}
		
		// add the final root index that will be removed later
		int root = tree.nodeList.size();
		// connect daughter and son node to temporary root node
		InsertRoot(tree, root);
		
		AssignThirdConnectedNode(tree);
		
		ParsimonyTree switchedTree0 = SwitchSubtrees(tree, internalNodeA, internalNodeB, 0);
		switchedTree0.PrintEdges(false, false, sb);
		sb.append(System.getProperty("line.separator"));
		ParsimonyTree switchedTree1 = SwitchSubtrees(tree, internalNodeA, internalNodeB, 1);
		switchedTree1.PrintEdges(false, false, sb);
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
	}
	
	/**
	 * This method assigns the parent internal node of leaf nodes.
	 * If the current node being checked is also an internal node, this method simply connects it to the third node
	 * which it connects to but haven't been registered yet as a connection.
	 * @param tree
	 */
	public void AssignThirdConnectedNode(ParsimonyTree tree) {
		int root = tree.nodeList.size()-1;
		for (int pIndex = tree.nodeList.size()-1; pIndex >= 0; pIndex--) {
			ParsimonyNode currNode = tree.nodeList.get(pIndex);
			if (currNode.isRoot()) {
				tree.nodeList.get(tree.tempRootInsPoint.nodeA).connection[2] = tree.nodeList.get(tree.tempRootInsPoint.nodeB);
				tree.nodeList.get(tree.tempRootInsPoint.nodeB).connection[2] = tree.nodeList.get(tree.tempRootInsPoint.nodeA);
			} else if (!currNode.isLeaf()) {
				ParsimonyNode left = currNode.connection[0];
				ParsimonyNode right = currNode.connection[1];
				left.connection[2] = currNode;
				right.connection[2] = currNode;
			}
		}
	}
	
	/**
	 * 
	 * @param tree
	 * @param internalNodeA
	 * @param internalNodeB
	 * @param switchIndex
	 */
	public ParsimonyTree SwitchSubtrees(ParsimonyTree tree, int nodeA, int nodeB, int switchIndex) {
		ParsimonyTree newTree = tree.CopyTree();
		ParsimonyNode internalNodeA = newTree.nodeList.get(nodeA);
		ParsimonyNode internalNodeB = newTree.nodeList.get(nodeB);
		
		int pivotA = -1; 
		int pivotB = -1;
		
		// assign the real position of the connected internal node from the connection array
		for (int i = 0; i < 3; i++) {
			if (internalNodeA.connection[i] != null && internalNodeA.connection[i].index == nodeB) {
				pivotA = i;
			}
			if (internalNodeB.connection[i] != null &&  internalNodeB.connection[i].index == nodeA) {
				pivotB = i;
			}
		}
		
		int xIndex = (pivotA+2)%3;
		int yIndex = (pivotB+2)%3;
		int zIndex = (pivotB+1)%3;
		
		ParsimonyNode x = internalNodeA.connection[xIndex];
		ParsimonyNode y = internalNodeB.connection[yIndex];
		ParsimonyNode z = internalNodeB.connection[zIndex];
		
		if (switchIndex == 0) {	//WX-YZ pattern to WY-XZ
			internalNodeB.connection[yIndex] = x;
			internalNodeA.connection[xIndex] = y;
			ReplaceParent(x, internalNodeA, internalNodeB);
			ReplaceParent(y, internalNodeB, internalNodeA);
		} else { //WX-YZ pattern to WZ-YX
			internalNodeB.connection[zIndex] = x;
			internalNodeA.connection[xIndex] = z;
			ReplaceParent(x, internalNodeA, internalNodeB);
			ReplaceParent(z, internalNodeB, internalNodeA);
		}
		return newTree;
	}

	private void ReplaceParent(ParsimonyNode childNode, ParsimonyNode parentToReplace, ParsimonyNode replaceBy) {
		ReplaceParent(childNode, parentToReplace, replaceBy, false);
	}
	
	/**
	 * Replaces one of the connections of a target node
	 * This is used for switching nodes/subtrees in the nearest neighbor tree algorithm
	 * @param childNode - the target node
	 * @param parentToReplace - the connection that will be replaced
	 * @param replaceBy - the new connection that will be created
	 * @param withRearrangement - if true, reassigned the new connection to the last connection node
	 */
	private void ReplaceParent(ParsimonyNode childNode, ParsimonyNode parentToReplace, ParsimonyNode replaceBy, boolean withRearrangement) {
		for (int i = 0; i < 3; i++) {
			if (childNode.connection[i] != null && childNode.connection[i].index == parentToReplace.index) {
				if (withRearrangement && i != 2) {
					ParsimonyNode originalThirdNode = childNode.connection[2];
					childNode.connection[2] = replaceBy;
					childNode.connection[i] = originalThirdNode;
					//childNode.connection[i] = replaceBy;
				} else {
					childNode.connection[i] = replaceBy;
				}
				break;
			}
		}
	}
	
	public void LargeParsimonyHeuristic(List<String> lines) {
		int score = Parsimony.INFINITY;
		ParsimonyTree tree = UnrootedSmallParsimony(lines);
		int newScore = tree.GetMinimumScore(false);
		ParsimonyTree newTree = tree;
		tree.PrintNodes();
		while (newScore < score) {
			score = newScore;
			tree = newTree;
			int neighborScore;
			ParsimonyTree neighborTree;
			System.out.println("while runs...");
			for (String key : tree.edgeList.keySet()) {
				ParsimonyEdge currEdge = tree.edgeList.get(key);
				System.out.println("edge for runs..." +currEdge.isInternal +" " +currEdge.nodeA +" " +currEdge.nodeB);
				if (currEdge.isInternal) {
					for (int switchIndex = 0; switchIndex <= 1; switchIndex++) {
						dnaLength = 0; //reset
						tree.PrintNodes();
						AssignThirdConnectedNode(tree); //TODO
						ParsimonyTree switchedTree = SwitchSubtrees(tree, currEdge.nodeA, currEdge.nodeB, switchIndex);
						System.out.println("switching trees " +switchIndex +" " +currEdge.nodeA +" " +currEdge.nodeB);
						switchedTree.PrintNodes();
						neighborTree = UnrootedSmallParsimony(switchedTree.TreeToAdjacencyList());
						neighborScore = neighborTree.GetMinimumScore(false);
						if (neighborScore < newScore) {
							newScore = neighborScore;
							newTree = neighborTree;
						}
					}
				}
			}
			newTree.PrintEdges(false, true, sb);
		}
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
	}
}

class TemporaryRootInsertionPoint {
	int nodeA;
	int nodeB;
	
	public TemporaryRootInsertionPoint(int nodeA, int nodeB) {
		this.nodeA = nodeA;
		this.nodeB = nodeB;
	}
}

class ParsimonyTree {
	Map<Integer, ParsimonyNode> nodeList;
	int leafCount; 
	Map<String, ParsimonyEdge> edgeList;
	TemporaryRootInsertionPoint tempRootInsPoint;
	
	public ParsimonyTree(int leafCount) {
		nodeList = new HashMap<Integer,ParsimonyNode>();
		edgeList = new HashMap<String, ParsimonyEdge>();
		this.leafCount = leafCount;
	}
	
	/**
	 * @return a new tree with the exact same set of nodes
	 * Do not copy the edges because this copy will be used to interchange nodes (via NeighboringTree algorithm).
	 * The edgeList should be reconstructed after the node interchange.
	 */
	public ParsimonyTree CopyTree() {
		ParsimonyTree newTree = new ParsimonyTree(leafCount);
		// do not use List.putAll function since it will also copy the references instead of creating a new node
		for (int i = 0; i < nodeList.size(); i++) {
			newTree.nodeList.put(i, nodeList.get(i).CopyNode());
		}
		// reassign the connection array for each node so that it points to the newly created nodes of newTree
		for (int i = 0; i < nodeList.size(); i++) {
			ParsimonyNode oldNode = nodeList.get(i);
			ParsimonyNode newNode = newTree.nodeList.get(i);
			for (int j = 0; j < oldNode.connection.length; j++) {
				if (oldNode.connection[j] != null) {
					newNode.connection[j] = newTree.nodeList.get(oldNode.connection[j].index);
				}
			}
		}
		return newTree;
	}
	
	/**
	 * adds a new node to the tree
	 * @param nodeIndex - the identifier of the node
	 * @param isLeaf - true if the node has only one edge connected to it
	 * @param child - the daughter (left) node if the parent node is also new; 
	 * 				- the son (right) node if parent node already exists
	 * 				- the parent node (or a fellow internal node for unrooted trees) if both children already exist
	 * @param label - the DNA string that represents this node
	 */
	public void AddNode(int nodeIndex, int child, String label) {
		if (!nodeList.containsKey(nodeIndex)) {
			nodeList.put(nodeIndex, new ParsimonyNode(nodeIndex, label));
		} else if (label.length() != 0) {
			nodeList.get(nodeIndex).label = label;
		}
		
		if (!nodeList.containsKey(child) && child != -1) {
			nodeList.put(child, new ParsimonyNode(child, ""));
		}

		ParsimonyNode node = nodeList.get(nodeIndex);
		int connIndex = (node.connection[0] == null) ? 0 : ((node.connection[1]==null) ? 1 : 2);
		if (node.IsNewConnection(child)) { 
			node.connection[connIndex] = nodeList.get(child);
		}
		
	}
	
	/**
	 * Input: An integer n followed by an adjacency list for a rooted binary tree with n leaves
	 * labeled by DNA strings.
     * Output: The minimum parsimony score of this tree, followed by the adjacency list of the
     * tree corresponding to labeling internal nodes by DNA strings in order to minimize the
     * parsimony score of the tree.
	 */
	public void SmallParsimony() {
		for (int i = 0; i < nodeList.size(); i++) {
			ParsimonyNode currNode = nodeList.get(i);
			if (currNode.isLeaf()) {
				currNode.tag = 1;
				for (int k = 0; k < Parsimony.ALPHABET.length; k++) {
					if (currNode.label.charAt(0) == Parsimony.ALPHABET[k]) {
						currNode.score[k] = 0;
					} else {
						currNode.score[k] = Parsimony.INFINITY;
					}
				}
			}
		}
		
		List<Integer> ripeNodes = GetRipeNodes();
		
		while (!ripeNodes.isEmpty()) {
			boolean isRoot = ripeNodes.size() == 1;
			for (int r = 0; r < ripeNodes.size(); r++) {
				ParsimonyNode ripeNode = nodeList.get(ripeNodes.get(r));
				ripeNode.tag = 1;
				ScorePair sp = ripeNode.GetEffectiveScore(isRoot);
				ripeNode.label = new String(sp.nucleotides);
			}
			ripeNodes = GetRipeNodes();
		}
		
		Set<Integer> listOfFixedLeaves = new HashSet<Integer>();
		for (int i = nodeList.size()-1; i >= leafCount; i--) {
			ParsimonyNode parent = nodeList.get(i);
			for (int j = 0; j < parent.connection.length; j++) {
				if (parent.connection[j] == null) {
					continue;
				}
				int childIndex = parent.connection[j].index;
				if (!listOfFixedLeaves.contains(childIndex)) {
					ChangeLabel(parent, nodeList.get(childIndex));
					listOfFixedLeaves.add(parent.connection[j].index);
				}				
			}
		}
	}
	
	/**
	 * If a child node has 2 minimum scores, it selects the nucleotide that will generate the minimum hamming 
	 * distance by copying the nucleotide from the parent. If the nucleotide in the parent is not among the 
	 * choices of nucleotides that generate the minimum score, the first choice is selected.
	 * @param parent
	 * @param child
	 */
	private void ChangeLabel(ParsimonyNode parent, ParsimonyNode child) {
		boolean isSet = false;
		if (child.label.length() > 1 && !child.isLeaf()) {
			char parentChar = parent.label.charAt(0);
			char[] choices = child.label.toCharArray();
			for (int c = 0; c < choices.length; c++) {
				if (choices[c] == parentChar) {
					child.label = parent.label.charAt(0) +"";
					isSet = true;
					break;
				}
			}
			if (!isSet) {
				child.label = child.label.substring(0,1);
			}
		}
	}
	
	/**
	 * @return a list of ripe nodes at the current iteration of the while loop in the parsimony algorithm
	 * Ripe node is defined as any node whose tag is 0 (unprocessed) but whose child nodes' tag is 1
	 */
	private List<Integer> GetRipeNodes() {
		List<Integer> ripeIndices = new ArrayList<Integer>();
		for (int i = 0; i < nodeList.size(); i++) {
			if (nodeList.get(i).isRipe()) {
				ripeIndices.add(i);
			}
		}
		return ripeIndices;
	}
	
	public void PrintNodes() {
		for (int i = 0; i < nodeList.size(); i++) {
			System.out.println(nodeList.get(i));
		}
	}
	
	/**
	 * gets the minimum score of the root among all nucleotide alphabet (A, C, T, G)
	 * This is computed by adding the hamming distances of all edges in the tree
	 * If the node is unrooted, the score of the edges that were originally connected to the root
	 * is counted twice so we need to subtract the length for one of the edges
	 */
	public int GetMinimumScore(boolean isRooted) {
		ComputeEdges(isRooted);
		Set<Integer> traversedNodes = new HashSet<Integer>();
		int score = nodeList.get(nodeList.size()-1).GetNodeScore(traversedNodes);
		if (!isRooted) {
			score -= nodeList.get(tempRootInsPoint.nodeA).edgeScore;
		}
		return score;
	}
	
	/**
	 * @param nodeA
	 * @param nodeB
	 */
	private String AddEdge(ParsimonyNode nodeA, ParsimonyNode nodeB) {
		String key = nodeA.index +BioinformaticsCommon.NODE_SEPARATOR +nodeB.index;
		if (!edgeList.containsKey(key)) {
			String nodeALabel = nodeA.label == "" ? Integer.toString(nodeA.index) : nodeA.label;
			String nodeBLabel = nodeB.label == "" ? Integer.toString(nodeB.index) : nodeB.label;
			ParsimonyEdge newEdge = new ParsimonyEdge(nodeALabel, nodeBLabel, nodeA.index, nodeB.index);
			newEdge.isInternal = !nodeA.isLeaf() && !nodeB.isLeaf();
			edgeList.put(key, newEdge);
			
		} else {
			//System.out.println("duplicate: " +key);
		}
		return key;
	}
	
	/**
	 * creates the list of edges by connecting all parent nodes with their daughter and child nodes
	 * Parent nodes are the nodes whose index is greater than or equal to the leaf count and whose daughter and
	 * son nodes are not null
	 * @param isRooted - if false, merge the edges between the temporary root node that was added and its children
	 */
	public void ComputeEdges(boolean isRooted) {
		for (int i = 0; i < nodeList.size(); i++) {
			ParsimonyNode parent = nodeList.get(i);
			if (parent.isLeaf()) {
				continue;
			}
			String key;
			if (!isRooted && i == nodeList.size() - 1) {	//if temporary root node of unrooted tree
				key = AddEdge(parent.connection[0], parent.connection[1]);
				parent.connection[0].edgeScore = edgeList.get(key).GetHammingDist();
				key = AddEdge(parent.connection[1], parent.connection[0]);
				parent.connection[1].edgeScore = edgeList.get(key).GetHammingDist();
				// remove all edges connected to the temporary root
				RemoveEdgeAndInverse(parent, parent.connection[0]);
				RemoveEdgeAndInverse(parent, parent.connection[1]);
			} else {
				AddEdge(parent, parent.connection[0]);
				key = AddEdge(parent.connection[0], parent);
				parent.connection[0].edgeScore = edgeList.get(key).GetHammingDist();
				AddEdge(parent, parent.connection[1]);
				key = AddEdge(parent.connection[1], parent);
				parent.connection[1].edgeScore = edgeList.get(key).GetHammingDist();
				AddEdge(parent, parent.connection[2]);
				key = AddEdge(parent.connection[2], parent);
				parent.connection[2].edgeScore = edgeList.get(key).GetHammingDist();
			}
		}
	}
	
	/**
	 * removes both the edges representing the connection between the two input nodes
	 * This is primarily used for removing the connection between the temporary root and its children
	 * @param nodeA
	 * @param nodeB
	 */
	private void RemoveEdgeAndInverse(ParsimonyNode nodeA, ParsimonyNode nodeB) {
		edgeList.remove(nodeA.index +BioinformaticsCommon.NODE_SEPARATOR +nodeB.index);
		edgeList.remove(nodeB.index +BioinformaticsCommon.NODE_SEPARATOR +nodeA.index);
	}
	
	/**
	 * prints all edges in the form a->b:c where a and b are parent and child nodes and c is their hamming distance
	 * Each pair of nodes will be represented as 2 edges (a->b and b->a)
	 * @param isRooted - whether the root of the tree is known or not
	 */
	public StringBuffer PrintEdges(boolean isRooted, boolean includeScore, StringBuffer sb) {
		int minParsimonyScore = GetMinimumScore(isRooted);
		if (includeScore) {
			sb.append(minParsimonyScore);
			sb.append(System.getProperty("line.separator"));
		}

		for (String key : edgeList.keySet()) {
			sb.append(edgeList.get(key));
			sb.append(System.getProperty("line.separator"));
		}
		return sb;
	}
	
	/**
	 * Convert back a tree to its adjacency list format.
	 * This will be used by the NeighborTree search algorithm that accepts an adjacency list as input parameter.
	 * @return the adjacency list that represents this tree (internal nodes shall be represented as indices
	 */
	public List<String> TreeToAdjacencyList() {
		List<String> lines = new ArrayList<String>();
		lines.add(Integer.toString(this.leafCount));
		ComputeEdges(false);
		for (String key : edgeList.keySet()) {
			ParsimonyEdge currEdge = edgeList.get(key);
			ParsimonyNode nodeA = nodeList.get(currEdge.nodeA);
			ParsimonyNode nodeB = nodeList.get(currEdge.nodeB);
			if (nodeA.isRoot() || nodeB.isRoot()) {
				continue;
			}
			String line = (nodeA.isLeaf() ? nodeA.label : nodeA.index) +BioinformaticsCommon.NODE_SEPARATOR +(nodeB.isLeaf() ? nodeB.label : nodeB.index);
			lines.add(line);
		}
		return lines;
	}
}

class ParsimonyEdge {
	String left;
	String right;
	int nodeA;
	int nodeB;
	int hammingDist;
	boolean isInternal;
	
	public ParsimonyEdge(String left, String right, int nodeA, int nodeB) {
		this.left = left;
		this.right = right;
		this.nodeA = nodeA;
		this.nodeB = nodeB;
	}
	
	public int GetHammingDist() {
		hammingDist = 0;
		for (int i = 0; i < left.length(); i++) {
			if (!left.substring(i, i+1).equals(right.substring(i, i+1))) {
				hammingDist++;	
			}
		}
		return hammingDist;
	}
	
	@Override
	public String toString() {
		return left +BioinformaticsCommon.NODE_SEPARATOR +right +(Parsimony.isScoreShown ? ":" +GetHammingDist() : "");
	}
}

class ParsimonyNode {
	int index;
	ParsimonyNode[] connection = new ParsimonyNode[3]; 
	int tag;
	String label;
	int[] score = new int[4];
	int edgeScore;
	
	public ParsimonyNode(int index, String label) {
		this.index = index;
		this.label = label;
	}
	
	public boolean IsNewConnection(int nodeIndex) {
		if ((connection[0] != null && connection[0].index == nodeIndex) 
				|| (connection[1] != null && connection[1].index == nodeIndex)
				|| (connection[2] != null && connection[2].index == nodeIndex))
			return false;
		return true;
	}
	
	/**
	 * @return a basic copy of the node; after all nodes are reconstructed in the CopyTree function of ParsimonyTree,
	 * reassign the connections using ParsimonyTree.CopyConnections(sourceNode, destNode). Simply setting the connection 
	 * of nodeCopy to the connection array of this node will just copy the references of this node and not the 
	 * new node created.
	 */
	public ParsimonyNode CopyNode() {
		ParsimonyNode nodeCopy = new ParsimonyNode(index, label);
		nodeCopy.tag = this.tag;
		nodeCopy.score = this.score;
		nodeCopy.edgeScore = this.edgeScore;
		return nodeCopy;
	}
	
	private int NonNullConnectionCount() {
		int count = 0;
		for (int i = 0; i < connection.length; i++) {
			if (connection[i] != null) {
				count++;
			}
		}
		return count;
	}

	public boolean isLeaf() {
		return NonNullConnectionCount() == 1;
	}
	
	public boolean isRoot() {
		return NonNullConnectionCount() == 2;
	}
	
	public boolean isRipe() {
		for (int i = 0; i < connection.length; i++) {
			if (connection[0] != null && !connection[0].isTagged()) {
				return false;
			}
		}
		return tag == 0;
	}
	
	/**
	 * if non-leaf, mark as tagged without checking for the value of tag
	 * if leaf, only mark as tagged if value of tag is 1
	 **/
	public boolean isTagged() {
		return (isLeaf() && tag == 1) || (!isLeaf());
	}
	
	public int GetNodeScore(Set<Integer> traversedNodes) {
		traversedNodes.add(index);
		if (isLeaf()) {
			return 0;
		} else {
			int runningScore = 0;
			for (int i = 0; i < connection.length; i++) {
				if (connection[i] != null && !traversedNodes.contains(connection[i].index)) {
					runningScore = runningScore + connection[i].edgeScore + connection[i].GetNodeScore(traversedNodes);
					traversedNodes.add(connection[i].index);
				}
			}
			return runningScore;
		}
	}
	
	public ScorePair GetEffectiveScore(boolean isRoot) {
		int[] daughterScores = connection[0].GetComponentScoreForParent();
		int[] sonScores = connection[1].GetComponentScoreForParent();
		int minScore = Parsimony.INFINITY;
		String selNucleotides = "";
		for (int i = 0; i < Parsimony.ALPHABET.length; i++) {
			score[i] = daughterScores[i] + sonScores[i];
			if (score[i] < minScore) {
				minScore = score[i];
				//selNucleotide = Parsimony.ALPHABET[i];
			}
			//System.out.print(score[i] +" ");
		}
		for (int i = 0; i < Parsimony.ALPHABET.length; i++) {
			if (score[i] == minScore) {
				selNucleotides += Parsimony.ALPHABET[i];
			}
		}
		/*
		 *  if node is the root node, ignore multiple values for the nucleotide at the current index
		 *  just select the first element
		 */
		if (isRoot) {
			selNucleotides = selNucleotides.substring(0, 1);
		}
		return new ScorePair(minScore, selNucleotides.toCharArray()); 
	}
	
	/**
	 * This method is used to compute the daughter and son component scores of a parent node
	 * @return the minimum score over all symbols for this node (adding a delta penalty of 1 for a mismatch)
	 */
	public int[] GetComponentScoreForParent() {
		int[] newScore = new int[4];
		int minScore = Parsimony.INFINITY;
		for (int i = 0; i < Parsimony.ALPHABET.length; i++) {
			if (score[i] < minScore) {
				minScore = score[i];
			}
		}
		for (int i = 0; i < Parsimony.ALPHABET.length; i++) {
			newScore[i] = minScore;
			if (score[i] > minScore) {
				newScore[i] += 1;
			}
		}
		return newScore;
	}
	
	
	@Override
	public String toString() {
		return index +" " +(label==""?"":label) +" " +(connection[0]==null?-1:connection[0].index) +" " +(connection[1]==null?-1:connection[1].index) +" " +(connection[2]==null?-1:connection[2].index);
	}
}

class ScorePair {
	int score;
	char[] nucleotides;
	
	public ScorePair(int score, char[] nucleotides) {
		this.score = score;
		this.nucleotides = nucleotides;
	}
}