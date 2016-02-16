import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
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
	
	public void UnrootedSmallParsimony(List<String> lines) {
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
			} else if (!mainTree.nodeList.containsKey(parentIndex) || mainTree.nodeList.get(parentIndex).connection[1] == null) {
				/* do not add the nodes who are already parents because the edge between these pair of nodes is
				 * where we will inject the temporary root node
				 */
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
		mainTree.AddNode(parentIndex, parentIndex - 1, "");
		mainTree.AddNode(parentIndex, parentIndex - 2, "");
		for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
			treeList.get(charIndex+1).AddNode(parentIndex, parentIndex - 1, "");
			treeList.get(charIndex+1).AddNode(parentIndex, parentIndex - 2, "");
		}
		
		for (int tIndex = 0; tIndex < treeList.size(); tIndex++) {
			ParsimonyTree tree = treeList.get(tIndex);
			for (int pIndex = tree.nodeList.size()-1; pIndex >= leafCount; pIndex--) {
				ParsimonyNode currNode = tree.nodeList.get(pIndex);
				ParsimonyNode left = currNode.connection[0];
				ParsimonyNode right = currNode.connection[1];
				left.connection[2] = currNode;
				right.connection[2] = currNode;
			}
		}
					
		//mainTree.PrintNodes();
		for (int k = 1; k <= dnaLength; k++) {
			treeList.get(k).SmallParsimony();
			for (int m = leafCount; m < mainTree.nodeList.size(); m++) {
				mainTree.nodeList.get(m).label += treeList.get(k).nodeList.get(m).label;
			}
		}
		mainTree.PrintEdges(false, true, sb);
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
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
		tree.AddNode(root, root - 1, "");
		tree.AddNode(root, root - 2, "");
		
		for (int pIndex = tree.nodeList.size()-1; pIndex >= 0; pIndex--) {
			ParsimonyNode currNode = tree.nodeList.get(pIndex);
			if (currNode.isRoot()) {
				tree.nodeList.get(root-1).connection[2] = tree.nodeList.get(root-2);
				tree.nodeList.get(root-2).connection[2] = tree.nodeList.get(root-1);
			} else if (!currNode.isLeaf()) {
				ParsimonyNode left = currNode.connection[0];
				ParsimonyNode right = currNode.connection[1];
				left.connection[2] = currNode;
				right.connection[2] = currNode;
			}
		}
		ParsimonyTree switchedTree0 = SwitchSubtrees(tree, internalNodeA, internalNodeB, 0);
		switchedTree0.PrintEdges(false, false, sb);
		sb.append(System.getProperty("line.separator"));
		ParsimonyTree switchedTree1 = SwitchSubtrees(tree, internalNodeA, internalNodeB, 1);
		switchedTree1.PrintEdges(false, false, sb);
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
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
			if (internalNodeA.connection[i].index == nodeB) {
				pivotA = i;
			}
			if (internalNodeB.connection[i].index == nodeA) {
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
		for (int i = 0; i < 3; i++) {
			if (childNode.connection[i] != null && childNode.connection[i].index == parentToReplace.index) {
				childNode.connection[i] = replaceBy;
				break;
			}
		}
	}
}

class ParsimonyTree {
	Map<Integer, ParsimonyNode> nodeList;
	int leafCount; 
	Map<String, ParsimonyEdge> edgeList;
	
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
		newTree.nodeList.putAll(nodeList);
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
		} 

		ParsimonyNode node = nodeList.get(nodeIndex);
		int connIndex = (node.connection[0] == null) ? 0 : ((node.connection[1]==null) ? 1 : 2);
		node.connection[connIndex] = nodeList.get(child);
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
		
		for (int i = nodeList.size()-1; i >= leafCount; i--) {
			ParsimonyNode parent = nodeList.get(i);
			ParsimonyNode daughter = parent.connection[0];
			ParsimonyNode son = parent.connection[1];
			ChangeLabel(parent, daughter);
			ChangeLabel(parent, son);
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
		int score = nodeList.get(nodeList.size()-1).GetNodeScore();
		if (!isRooted) {
			score -= nodeList.get(nodeList.size()-2).edgeScore;
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
			edgeList.put(key, new ParsimonyEdge(nodeALabel, nodeBLabel));
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
	 * prints all edges in the form a->b:c where a and b are parent and child nodes and c is their hamming distance
	 * Each pair of nodes will be represented as 2 edges (a->b and b->a)
	 * @param isRooted - whether the root of the tree is known or not
	 */
	public StringBuffer PrintEdges(boolean isRooted, boolean includeScore, StringBuffer sb) {
		ComputeEdges(isRooted);
		if (includeScore) {
			sb.append(GetMinimumScore(isRooted));
			sb.append(System.getProperty("line.separator"));
		}
		// commented code is for sorting by left node (the a in a -> b)
		// this is not needed but makes the list easier to check when there are a lot of nodes
		/*Comparator<ParsimonyEdge> byNodeA = new Comparator<ParsimonyEdge>() {
	        @Override
	        public int compare(ParsimonyEdge edge1, ParsimonyEdge edge2) {
	        	return edge1.left.compareTo(edge2.left);
	        }
	    };
		Collections.sort(edgeList, byNodeA);*/
		for (String key : edgeList.keySet()) {
			sb.append(edgeList.get(key));
			sb.append(System.getProperty("line.separator"));
		}
		return sb;
	}
}

class ParsimonyEdge {
	String left;
	String right;
	int hammingDist;
	
	public ParsimonyEdge(String left, String right) {
		this.left = left;
		this.right = right;
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

class LargeParsimonyNode {
	
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

	public boolean isLeaf() {
		return connection[1] == null;
	}
	
	public boolean isRoot() {
		return connection[2] == null;
	}
	
	public boolean isRipe() {
		return tag == 0 && connection[0].tag == 1 && connection[1].tag == 1;
	}
	
	public int GetNodeScore() {
		if (isLeaf()) {
			return 0;
		} else {
			return connection[0].edgeScore + connection[1].edgeScore + connection[0].GetNodeScore() + connection[1].GetNodeScore();
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