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
	
	/**
	 * performs the small parsimony algorithm on each of the character indices of a rooted tree
	 * (i.e. if there are 8 characters on the dna string, this algorithm is performed 8 times and then the
	 * all the results are combined together to form the dna string for each internal node (and root) that 
	 * produces the minimum parsimony score for the input adjacency list)
	 * @param lines - a number on the first line representing the number of leaves in the tree followed by
	 * 					an adjacency list representing the predefined rooted tree structure
	 */
	public void RootedParsimony(List<String> lines) {
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
				mainTree.AddNode(i-1, true, -1, childLabel);
				mainTree.AddNode(parentIndex, false, i-1, "");
				
				for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
					treeList.get(charIndex+1).AddNode(i-1, true, -1, childLabel.substring(charIndex, charIndex+1));
					treeList.get(charIndex+1).AddNode(parentIndex, false, i-1, "");
				}
			} else {
				int childIndex = Integer.parseInt(st.nextToken());
				mainTree.AddNode(parentIndex, false, childIndex, "");
				for (int charIndex = 0; charIndex < dnaLength; charIndex++) {
					treeList.get(charIndex+1).AddNode(parentIndex, false, childIndex, "");
				}
			}
		}
		
		for (int k = 1; k <= dnaLength; k++) {
			treeList.get(k).SmallParsimony();
			for (int m = leafCount; m < mainTree.nodeList.size(); m++) {
				mainTree.nodeList.get(m).leafLabel += treeList.get(k).nodeList.get(m).leafLabel;
			}
		}
		mainTree.PrintEdges();
	}
	
	
	
	

}

class ParsimonyTree {
	Map<Integer, ParsimonyNode> nodeList;
	int leafCount; 
	List<ParsimonyEdge> edgeList;
	
	public ParsimonyTree(int leafCount) {
		nodeList = new HashMap<Integer,ParsimonyNode>();
		edgeList = new ArrayList<ParsimonyEdge>();
		this.leafCount = leafCount;
	}
	
	/**
	 * adds a new node to the tree
	 * @param nodeIndex - the identifier of the node
	 * @param isLeaf - true if the node has only one edge connected to it
	 * @param child - the daughter (left) node if the parent node is also new; 
	 * 				- the son (right) node if parent node already exists
	 * @param leafLabel - the DNA string that represents this node
	 */
	public void AddNode(int nodeIndex, boolean isLeaf, int child, String leafLabel) {
		if (!nodeList.containsKey(nodeIndex)) {
			nodeList.put(nodeIndex, new ParsimonyNode(nodeIndex, isLeaf, nodeList.get(child), null, leafLabel));
		} else {
			ParsimonyNode node = nodeList.get(nodeIndex);
			node.son = nodeList.get(child);
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
			if (currNode.isLeaf) {
				currNode.tag = 1;
				for (int k = 0; k < Parsimony.ALPHABET.length; k++) {
					if (currNode.leafLabel.charAt(0) == Parsimony.ALPHABET[k]) {
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
				ripeNode.leafLabel = new String(sp.nucleotides);
			}
			ripeNodes = GetRipeNodes();
		}
		
		for (int i = nodeList.size()-1; i >= leafCount; i--) {
			ParsimonyNode parent = nodeList.get(i);
			ParsimonyNode daughter = parent.daughter;
			ParsimonyNode son = parent.son;
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
		if (child.leafLabel.length() > 1 && !child.isLeaf) {
			char parentChar = parent.leafLabel.charAt(0);
			char[] choices = child.leafLabel.toCharArray();
			for (int c = 0; c < choices.length; c++) {
				if (choices[c] == parentChar) {
					child.leafLabel = parent.leafLabel.charAt(0) +"";
					isSet = true;
					break;
				}
			}
			if (!isSet) {
				child.leafLabel = child.leafLabel.substring(0,1);
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
	 */
	public int GetMinimumScore() {
		return nodeList.get(nodeList.size()-1).GetNodeScore();
	}
	
	/**
	 * creates the list of edges by connecting all parent nodes with their daughter and child nodes
	 * Parent nodes are the nodes whose index is greater than or equal to the leaf count and whose daughter and
	 * son nodes are not null
	 */
	public void ComputeEdges() {
		for (int i = leafCount; i < nodeList.size(); i++) {
			ParsimonyNode parent = nodeList.get(i);
			edgeList.add(new ParsimonyEdge(parent.leafLabel, parent.daughter.leafLabel));
			edgeList.add(new ParsimonyEdge(parent.daughter.leafLabel, parent.leafLabel));
			parent.daughter.edgeScore = edgeList.get(edgeList.size()-1).GetHammingDist();
			edgeList.add(new ParsimonyEdge(parent.leafLabel, parent.son.leafLabel));
			edgeList.add(new ParsimonyEdge(parent.son.leafLabel, parent.leafLabel));
			parent.son.edgeScore = edgeList.get(edgeList.size()-1).GetHammingDist();
		}
	}
	
	/**
	 * prints all edges in the form a->b:c where a and b are parent and child nodes and c is their hamming distance
	 * Each pair of nodes will be represented as 2 edges (a->b and b->a)
	 */
	public void PrintEdges() {
		ComputeEdges();
		StringBuffer sb = new StringBuffer();
		sb.append(GetMinimumScore());
		sb.append(System.getProperty("line.separator"));
		// commented code is for sorting by left node (the a in a -> b)
		// this is not needed but makes the list easier to check when there are a lot of nodes
		/*Comparator<ParsimonyEdge> byNodeA = new Comparator<ParsimonyEdge>() {
	        @Override
	        public int compare(ParsimonyEdge edge1, ParsimonyEdge edge2) {
	        	return edge1.left.compareTo(edge2.left);
	        }
	    };
		Collections.sort(edgeList, byNodeA);*/
		for (int i = 0; i < edgeList.size(); i++) {
			sb.append(edgeList.get(i));
			sb.append(System.getProperty("line.separator"));
		}
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
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
		return left +BioinformaticsCommon.NODE_SEPARATOR +right +":" +GetHammingDist();
	}
}

class ParsimonyNode {
	boolean isLeaf;
	int index;
	ParsimonyNode daughter;
	ParsimonyNode son;
	int tag;
	String leafLabel;
	int[] score = new int[4];
	int edgeScore;
	
	public ParsimonyNode(int index, boolean isLeaf, ParsimonyNode daughter, ParsimonyNode son, String leafLabel) {
		this.index = index;
		this.isLeaf = isLeaf;
		this.daughter = daughter;
		this.son = son;
		this.leafLabel = leafLabel;
	}
	
	public boolean isRipe() {
		return tag == 0 && daughter.tag == 1 && son.tag == 1;
	}
	
	public int GetNodeScore() {
		if (isLeaf) {
			return 0;
		} else {
			return daughter.edgeScore + son.edgeScore + daughter.GetNodeScore() + son.GetNodeScore();
		}
	}
	
	public ScorePair GetEffectiveScore(boolean isRoot) {
		int[] daughterScores = daughter.GetComponentScoreForParent();
		int[] sonScores = son.GetComponentScoreForParent();
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
		return index +" " +(leafLabel==""?"":leafLabel) +" " +(daughter==null?-1:daughter.index) +" " +(son==null?-1:son.index);
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