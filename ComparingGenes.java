import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

public class ComparingGenes {

	public static void main(String[] args) {
		List<String> params = BioinformaticsCommon.readFile();

		// test for multiple sequence alignment
		String v = params.get(0);
		String w = params.get(1);
		String u = params.get(2);
		MultipleLongestCommonSubsequence(v, w, u);
	}
	
	/**
	 * @param money - the amount of money that is being translated to coins
	 * @param coins - an array containing the denominations
	 */
	private static void DPChange(int money, List<Integer> coins) {
		int[] minCoins = new int[money+1];
		for(int i = 1; i <= money; i++) {
			minCoins[i] = Integer.MAX_VALUE;
			for (int j = 0; j < coins.size(); j++) {
				int coinValue = coins.get(j);
				if (i >= coinValue) {
					if (minCoins[i-coinValue] + 1 < minCoins[i]) {
						minCoins[i] = minCoins[i-coinValue] + 1;
					}
				}
			}
		}
		System.out.println(minCoins[money]);
	}
	
	/**
	 * Input: Integers n and m, followed by an n × (m + 1) matrix Down and an (n + 1) × m matrix Right.
     * The two matrices are separated by the - symbol.
     * Output: The length of a longest path from source (0, 0) to sink (n, m) in the n × m rectangular grid
     * whose edges are defined by the matrices Down and Right.
	 */
	private static void ManhattanTouristLongestPath(int n, int m, int[][] down, int[][] right) {
		int s[][] = new int[n+1][m+1];	//the nodes in the manhattan tourist matrix
		s[0][0] = 0;
		for (int i = 1; i <= n; i++) {
			s[i][0] = s[i-1][0] + down[i-1][0];
		}
		for (int j = 1; j <= m; j++) {
			s[0][j] = s[0][j-1] + right[0][j-1];
		}
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= m; j++) {
				s[i][j] = Math.max(s[i-1][j] + down[i-1][j], s[i][j-1] + right[i][j-1]);
			}
		}
		System.out.println(s[n][m]);
	}
	
	private static final int DELETION = 0;
	private static final int INSERTION = 1;
	private static final int MATCH = 2;
	private static final int JUMP_PATH = 3;
	
	/**
	 * 	 Input: Two strings s and t.
	 *   Output: A longest common subsequence of s and t. (Note: more than one solution may exist,
	 *    in which case you may output any one.)
	 */     
	private static void LCSBacktrack(String v, String w) {
		int rows = v.length() + 1;
		int cols = w.length() + 1;
		int s[][] = new int[rows][cols];
		int backtrack[][] = new int[rows][cols];
		
		for (int i = 0; i < rows; i++) {
			s[i][0] = 0;
		}
		for (int j = 0; j < cols; j++) {
			s[0][j] = 0;
		}
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				int ins = s[i][j-1];
				int del = s[i-1][j];
				boolean isMatching = v.charAt(i-1) == w.charAt(j-1);
				int match = s[i-1][j-1] + (isMatching ? 1 : 0);
				s[i][j] = Math.max(Math.max(ins, del), match);
				if (s[i][j] == del) {
					backtrack[i][j] = DELETION;
				} else if (s[i][j] == ins) {
					backtrack[i][j] = INSERTION;
				} else if (s[i][j] == s[i-1][j-1] + 1 && isMatching) {
					backtrack[i][j] = MATCH;
				}
			}
		}
		StringBuffer sb = new StringBuffer();
		OutputLCS(backtrack, v, rows-1, cols-1, sb);
		System.out.println(sb.toString());
		BioinformaticsCommon.WriteOutputToFile(sb.toString());
	}
	
	/**
	 * Utility method for "backtracking" through the longest path and getting the matching nucleotides along the way
	 * @param backtrack - the matrix containing "DELETION", "INSERTION" or "MATCH" for each element. 
	 * 						The content is determined by checking which method was used to get the maximum value
	 * 						for the specific node (with priority as in the order described above)
	 * @param v - the first dna string in the 2 dna strings being matched
	 * @param i - the row count of the current node being backtracked
	 * @param j - the column count of the current node being backtracked
	 * @param sb - a buffer that stores the nucleotide matches that are identified along the way
	 */
	private static void OutputLCS(int[][] backtrack, String v, int i, int j, StringBuffer sb) {
		if (i == 0 || j == 0)
			return;
		if (backtrack[i][j] == DELETION) {
			OutputLCS(backtrack, v, i-1, j, sb);
		} else if (backtrack[i][j] == INSERTION) {
			OutputLCS(backtrack, v, i, j-1, sb);
		} else if (backtrack[i][j] == MATCH) {
			OutputLCS(backtrack, v, i-1, j-1, sb);
			sb.append(v.charAt(i-1));
		}
	}
	
	private static void OutputLCSBasedOnScoringMatrix(int[][] backtrack, String v, String w, int i, int j, StringBuffer sbv, StringBuffer sbw) {
		OutputLCSBasedOnScoringMatrix(backtrack, v, w, i, j, sbv, sbw, AlignmentType.GLOBAL);
	}
	
	private static void OutputLCSBasedOnScoringMatrix(int[][] backtrack, String v, String w, int i, int j, StringBuffer sbv, StringBuffer sbw, AlignmentType alignType) {
		if (i == 0 || j == 0) {
			if (alignType.equals(AlignmentType.GLOBAL)) {
				if (i < j) {
					sbv.insert(0, String.format("%" +(j) +"s", "").replace(' ', '-'));
					sbw.insert(0, w.substring(i, j));
				} else if (i > j) {
					sbv.insert(0, v.substring(j, i));
					sbw.insert(0, String.format("%" +(i) +"s", "").replace(' ', '-'));
				}
			}
			return;
		}
		//System.out.println(i +" " +j +" " +backtrack[i][j]);
		if (backtrack[i][j] == DELETION) {
			OutputLCSBasedOnScoringMatrix(backtrack, v, w, i-1, j, sbv, sbw, alignType);
			sbv.append(v.charAt(i-1));
			sbw.append("-");
		} else if (backtrack[i][j] == INSERTION) {
			OutputLCSBasedOnScoringMatrix(backtrack, v, w, i, j-1, sbv, sbw, alignType);
			sbv.append("-");
			sbw.append(w.charAt(j-1));
		} else if (backtrack[i][j] == MATCH) {
			OutputLCSBasedOnScoringMatrix(backtrack, v, w, i-1, j-1, sbv, sbw, alignType);
			sbv.append(v.charAt(i-1));
			sbw.append(w.charAt(j-1));
		} else { //ZERO_PATH
			//do nothing;
		}
	}
	
	/**
	 *  Input: An integer representing the source node of a graph, followed by an integer representing the
     *  sink node of the graph, followed by a list of edges in the graph. The edge notation 0->1:7 indicates
     *  that an edge connects node 0 to node 1 with weight 7. 
     *  Output: The length of a longest path in the graph, followed by a longest path. (If multiple longest paths exist, you may return any one.)
	 */
	private static void LongestPathDAG(int source, int sink, List<String> edges) {
		List<Node> nodes = new ArrayList<Node>();
		for (int i = 0; i < edges.size(); i++) {
			StringTokenizer st = new StringTokenizer(edges.get(i), "->");
			int predecessor = Integer.parseInt(st.nextToken());
			StringTokenizer st2 = new StringTokenizer(st.nextToken(), ":");
			int currNodeIndex = Integer.parseInt(st2.nextToken());
			int currNodeWeight = Integer.parseInt(st2.nextToken());
			
			if (!nodes.contains(new Node(currNodeIndex))) 
				nodes.add(new Node(currNodeIndex));
			if (!nodes.contains(new Node(predecessor))) 
				nodes.add(new Node(predecessor));
			
			nodes.get(nodes.indexOf(new Node(currNodeIndex))).addPredecessor(nodes.get(nodes.indexOf(new Node(predecessor))), currNodeWeight);
		}
		Collections.sort(nodes);
		for (int m = 0; m < nodes.size(); m++) {
			nodes.get(m).computeMaxValue(source);
		}
		List<Integer> backtrackList = new ArrayList<Integer>();
		backtrackList.add(sink);
		int backtrackPredecessor = nodes.get(nodes.indexOf(new Node(sink))).backtrackPredecessor.get(0);
		while (backtrackPredecessor != source) {
			backtrackList.add(backtrackPredecessor);
			backtrackPredecessor = nodes.get(nodes.indexOf(new Node(backtrackPredecessor))).backtrackPredecessor.get(0);
		}
		backtrackList.add(source);
		String output = BioinformaticsCommon.joinListReverse(backtrackList, "->");
		System.out.println(nodes.get(nodes.indexOf(new Node(sink))).value);
		System.out.println(output);
	}
	
	/**
	 * 
	 * @param v - the first dna string to be matched
	 * @param w - the second dna string to be matched
	 * @param scoreMatrix - the protein scoring matrix (e.g. BLOSUM62)
	 * @param mismatchPenalty - this is only used when score matrix is null (i.e. penalties are fixed)
	 * @param indelPenalty - the penalty for insertions and deletions
	 * @param alignType - value is either global, local, fitting, overlap
	 * @return the Alignment object containing the maximum score, aligned first dna string and aligned second dna string
	 * Input: Two protein strings written in the single-letter amino acid alphabet.
	 * Output: 
	 * 1.) If solving for global alignment: The maximum alignment score of these strings followed by an alignment achieving this
     * maximum score.
     * 2.) If solving for local alignment: The maximum score of a local alignment of the strings, followed by a local alignment of these
     * strings achieving the maximum score.
	 */
	private static Alignment LCSWithScoringMatrix(String v, String w, int[][] scoreMatrix, int mismatchPenalty, int indelPenalty, AlignmentType alignType) {
		int rows = v.length() + 1;
		int cols = w.length() + 1;
		int s[][] = new int[rows][cols];
		int backtrack[][] = new int[rows][cols];
		List<String> aminoList = BioinformaticsCommon.AMINO_ACID_LIST_FOR_SCORING;
		int maxLocalScore = 0;
		int maxLocalIndexI = rows - 1;
		int maxLocalIndexJ = cols - 1;
		
		for (int i = 0; i < rows; i++) {
			s[i][0] = 0 - (i*indelPenalty);
		}
		for (int j = 0; j < cols; j++) {
			s[0][j] = 0 - (j*indelPenalty);
		}
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				int ins = s[i][j-1] - indelPenalty;
				int del = s[i-1][j] - indelPenalty;
				String vStr = v.charAt(i-1) +"";
				String wStr = w.charAt(j-1) +"";
				int charScore = (vStr.equals(wStr) ? 1 : (-1*mismatchPenalty));
				if (scoreMatrix != null) {
					charScore = scoreMatrix[aminoList.indexOf(vStr)][aminoList.indexOf(wStr)];
				}
				int match = s[i-1][j-1] + charScore;
				if (alignType.equals(AlignmentType.GLOBAL)) {
					s[i][j] = Math.max(Math.max(ins, del), match);
				} else if (alignType.equals(AlignmentType.LOCAL)) {
					s[i][j] = Math.max(Math.max(Math.max(ins, del), match), 0);
					if (s[i][j] > maxLocalScore) {
						maxLocalIndexI = i;
						maxLocalIndexJ = j;
						maxLocalScore = s[i][j];
					}		
				} else if (alignType.equals(AlignmentType.FITTING)) {
					s[i][j] = Math.max(Math.max(ins, del), match);
					if (j == cols-1) {
						if (s[i][j] > maxLocalScore) {
							maxLocalIndexI = i;
							maxLocalIndexJ = cols-1;
							maxLocalScore = s[i][j];
						}
					}
				} else if (alignType.equals(AlignmentType.OVERLAP)) {
					s[i][j] = Math.max(Math.max(ins, del), match);
					if (i == rows-1) {
						if (s[i][j] >= maxLocalScore) {
							maxLocalIndexI = rows-1;
							maxLocalIndexJ = j;
							maxLocalScore = s[i][j];
						}
					}
				} 
				
				if (s[i][j] == del) {
					backtrack[i][j] = DELETION;
				} else if (s[i][j] == ins) {
					backtrack[i][j] = INSERTION;
				} else if (s[i][j] == match) {
					backtrack[i][j] = MATCH;
				} else { 
					backtrack[i][j] = JUMP_PATH;
				}
			}
		}
		
		StringBuffer sbv = new StringBuffer();
		StringBuffer sbw = new StringBuffer();
		OutputLCSBasedOnScoringMatrix(backtrack, v, w, maxLocalIndexI, maxLocalIndexJ, sbv, sbw, alignType);
		Alignment alignment = new Alignment(s[maxLocalIndexI][maxLocalIndexJ], sbv.toString(), sbw.toString());
		return alignment;
	}
	
	private static void LCSWithScoringMatrix(String v, String w, int[][] scoreMatrix, int mismatchPenalty, int indelPenalty) {
		LCSWithScoringMatrix(v, w, scoreMatrix, mismatchPenalty, indelPenalty, AlignmentType.GLOBAL);
	}
	
	/**
	 * Input: Two strings.
     * Output: The edit distance between these strings.
     *  		Edit distance is defined as the number of insertion, deletions and replacements needed to convert 
     *  		one string to another.
	 */
	int editDistance = 0;
	private static void FindEditDistance(String v, String w) {
		int[][] editDist = new int[v.length()][w.length()];
		for (int i = 0; i < v.length(); i++) {
			for (int j = 0; j < w.length(); j++) {
				if (i == 0) {
					editDist[i][j] = j;
				} else if (j == 0) {
					editDist[i][j] = i;
				} else if (v.charAt(i-1) == w.charAt(j-1)) {
					editDist[i][j] = editDist[i-1][j-1];
				} else {
					editDist[i][j] = 1 + Math.min(Math.min(editDist[i-1][j], editDist[i][j-1]), editDist[i-1][j-1]);
				}
			}
		}
		System.out.println(editDist[v.length()-1][w.length()-1]);
	}

	/**
	 * 
	 * @param v - the first dna/amino acid string to match
	 * @param w - the second dna/amino acid string to match
	 * @param scoreMatrix - the protein scoring matrix (e.g. BLOSUM62)
	 * @param sigma - the gap opening penalty (score deduction for the first indel in a series of indels)
	 * @param epsilon - the gap extension penalty (score deduction for succeeding continuous indels)
	 * Input: Two amino acid strings v and w (each of length at most 100).
     * Output: The maximum alignment score between v and w, followed by an alignment of v and w
     * achieving this maximum score. Use the BLOSUM62 scoring matrix, a gap opening penalty of 11, and
     * a gap extension penalty of 1.
     * Sample Input:
     * PRTEINS
     * PRTWPSEIN
	 * Sample Output:
     * 8
     * PRT---EINS
     * PRTWPSEIN-
	 */
	private static Alignment AlignmentWithAffineGapPenalties(String v, String w, int[][] scoreMatrix, int sigma, int epsilon) {
		int rows = v.length() + 1;
		int cols = w.length() + 1;
		
		//create three manhattan matrices
		int lower[][] = new int[rows][cols];
		int middle[][] = new int[rows][cols];
		int upper[][] = new int[rows][cols];
		
		int backtrack[][] = new int[rows][cols];
		List<String> aminoList = BioinformaticsCommon.AMINO_ACID_LIST_FOR_SCORING;
		int maxLocalIndexI = rows - 1;
		int maxLocalIndexJ = cols - 1;
		
		for (int i = 0; i < rows; i++) {
			lower[i][0] = 0;
			middle[i][0] = 0;
			upper[i][0] = 0;
		}
		for (int j = 0; j < cols; j++) {
			lower[0][j] = 0;
			middle[0][j] = 0;
			upper[0][j] = 0;
		}
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				String vStr = v.charAt(i-1) +"";
				String wStr = w.charAt(j-1) +"";
				int charMatchingScore = scoreMatrix[aminoList.indexOf(vStr)][aminoList.indexOf(wStr)];
				
				lower[i][j] = Math.max(lower[i-1][j] - epsilon, middle[i-1][j] - sigma);
				upper[i][j] = Math.max(upper[i][j-1] - epsilon, middle[i][j-1] - sigma);
				middle[i][j] = Math.max(middle[i-1][j-1] + charMatchingScore, Math.max(lower[i][j], upper[i][j]));
				
				
				if (middle[i][j] == lower[i][j]) {
					backtrack[i][j] = DELETION;
				} else if (middle[i][j] == upper[i][j]) {
					backtrack[i][j] = INSERTION;
				} else {
					backtrack[i][j] = MATCH;
				}
			}
		}
		
		StringBuffer sbv = new StringBuffer();
		StringBuffer sbw = new StringBuffer();
		OutputLCSBasedOnScoringMatrix(backtrack, v, w, maxLocalIndexI, maxLocalIndexJ, sbv, sbw, AlignmentType.GLOBAL);
		Alignment alignment = new Alignment(middle[maxLocalIndexI][maxLocalIndexJ], sbv.toString(), sbw.toString());
		return alignment;
	}
	
	
	/**
	 * 
	 * @param v - the first dna/amino acid string to match
	 * @param w - the second dna/amino acid string to match
	 * @param scoreMatrix - the protein scoring matrix (e.g. BLOSUM62)
	 * @param indelPenalty - a linear penalty for insertion and deletions
	 */
	private static MiddleElement MiddleEdgeInLinearSpace(String v, String w, int[][] scoreMatrix, int indelPenalty) {
		int rows = v.length() + 1;
		int cols = w.length() + 1;
		int middleCol = (cols-1)/2;
		int[] middleArray = new int[rows];
		int[] nextToMiddleArray = new int[rows];
		int bestMiddleScore = 0;
		int bestMiddleScoreRow = 0;
		
		int s[][] = new int[rows][middleCol+2];
		// contains the max score for the next node from the middle node (order: delete, match, insert)
		int[] nextEdgeScore = new int[3];
		
		List<String> aminoList = BioinformaticsCommon.AMINO_ACID_LIST_FOR_SCORING;
		
		for (int i = 0; i < rows; i++) {
			s[i][0] = 0 - (i * indelPenalty);
		}
		for (int j = 0; j <= middleCol+1; j++) {
			s[0][j] = 0 - (j * indelPenalty);
			if (j == middleCol) {
				middleArray[0] = s[0][j];
			} else if (j == middleCol+1) {
				nextToMiddleArray[0] = s[0][j];
			}
		}
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j <= middleCol+1; j++) {
				String vStr = v.charAt(i-1) +"";
				String wStr = w.charAt(j-1) +"";
				int charMatchingScore = vStr.equals(wStr) ? 1 : -1;
				if (scoreMatrix != null)
					charMatchingScore = scoreMatrix[aminoList.indexOf(vStr)][aminoList.indexOf(wStr)];
				s[i][j] = Math.max(s[i-1][j-1] + charMatchingScore, Math.max(s[i-1][j] - indelPenalty, s[i][j-1] - indelPenalty));
				if (j == middleCol) {
					middleArray[i] = s[i][j];
				} else if (j == middleCol + 1) {
					nextToMiddleArray[i] = s[i][j];
				}

			}
		}
		
		String vRev = reverseString(v);
		String wRev = reverseString(w);
		for (int i = 0; i < rows; i++) {
			s[i][0] = 0 - (i*indelPenalty);
		}
		for (int j = 0; j < middleCol; j++) {
			s[0][j] = 0 - (j*indelPenalty);
		}
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j <= middleCol; j++) {
				String vStr = vRev.charAt(i-1) +"";
				String wStr = wRev.charAt(j-1) +"";
				int charMatchingScore = vStr.equals(wStr) ? 1 : 0;
				if (scoreMatrix != null)
					charMatchingScore = scoreMatrix[aminoList.indexOf(vStr)][aminoList.indexOf(wStr)];
				s[i][j] = Math.max(s[i-1][j-1] + charMatchingScore, Math.max(s[i-1][j] - indelPenalty, s[i][j-1] - indelPenalty));
				if (j == middleCol) {
					//middleArray[rows-1-i] += s[i][j];
					nextToMiddleArray[rows-1-i] += s[i][j-1];
					if (middleArray[rows-1-i] >= bestMiddleScore) {
						bestMiddleScore = middleArray[rows-1-i];
						bestMiddleScoreRow = rows-1-i;
					}
				}

			}
		}
		
		int middleEdge;
		
		System.out.print("(" +bestMiddleScoreRow +", " +middleCol +") ");
		nextEdgeScore[0] = middleArray[bestMiddleScoreRow+1];
		nextEdgeScore[1] = nextToMiddleArray[bestMiddleScoreRow+1];
		nextEdgeScore[2] = nextToMiddleArray[bestMiddleScoreRow];
		int maxNextEdgeScore = Math.max(nextEdgeScore[0], Math.max(nextEdgeScore[1], nextEdgeScore[2]));
		if (maxNextEdgeScore == nextEdgeScore[1]) {
			System.out.print("(" +(bestMiddleScoreRow+1) +", " +(middleCol+1) +")");
			middleEdge = MATCH;
		} else if (maxNextEdgeScore == nextEdgeScore[0]) {
			System.out.print("(" +(bestMiddleScoreRow+1) +", " +middleCol +")");
			middleEdge = DELETION;
		} else {
			System.out.print("(" +bestMiddleScoreRow +", " +(middleCol+1) +")");
			middleEdge = INSERTION;
		}
		
		return new MiddleElement(middleCol, bestMiddleScoreRow, middleEdge);
		
	}
	
	/**
	 * @param s - the string to be reversed
	 * @return - the reversed string
	 */
	private static String reverseString(String s) {
		StringBuilder sb = new StringBuilder(s);
		return sb.reverse().toString();
	}
	
	/**
	 * used for aligning three strings (global alignment is used)
	 * @param v - the first string to be aligned
	 * @param w - the second string to be aligned
	 * @param u - the third string to be aligned
	 */
	private static void MultipleLongestCommonSubsequence(String v, String w, String u) {
		int rows = v.length() + 1;
		int cols = w.length() + 1;
		int breadth = u.length() + 1;
		int s[][][] = new int[rows][cols][breadth];
		int backtrack[][][] = new int[rows][cols][breadth];
		
		for (int i = 0; i > rows; i++) {
			s[i][0][0] = 0;
		}
		for (int j = 0; j < cols; j++) {
			s[0][j][0] = 0;
		}
		for (int k = 0; k < breadth; k++) {
			s[0][0][k] = 0;
		}
		
		for (int i = 1; i < rows; i++) {
			for (int j = 1; j < cols; j++) {
				for (int k = 1; k < breadth; k++) {
					String vSub = v.substring(i-1, i);
					String wSub = w.substring(j-1, j);
					String uSub = u.substring(k-1, k);
					
					int matchScore = vSub.equals(wSub) && vSub.equals(uSub) ? 1 : 0;
					
					int[] dims = new int[7];
					dims[0] = s[i-1][j][k];
					dims[1] = s[i][j-1][k];
					dims[2] = s[i][j][k-1];
					dims[3] = s[i-1][j-1][k];
					dims[4] = s[i-1][j][k-1];
					dims[5] = s[i][j-1][k-1];
					dims[6] = s[i-1][j-1][k-1] + matchScore;
					
					int[] copy = new int[7];
					System.arraycopy(dims, 0, copy, 0, dims.length);
					Arrays.sort(copy);
					s[i][j][k] = copy[6];
					
					for (int m = 0; m < 7; m++) {
						if (copy[6] == dims[m]) {
							backtrack[i][j][k] = m;
						}
					}
				}
			}
		}
		
		System.out.println(s[rows-1][cols-1][breadth-1]);
		StringBuffer sbv = new StringBuffer();
		StringBuffer sbw = new StringBuffer();
		StringBuffer sbu = new StringBuffer();
		ThreeDimensionalBacktrack(backtrack, rows-1, cols-1, breadth-1, sbv, sbw, sbu, v, w, u);
		System.out.println(sbv.reverse());
		System.out.println(sbw.reverse());
		System.out.println(sbu.reverse());
	}
	
	/**
	 * constructs the alignment strings of three strings by going from sink to source
	 * @param backtrack - the backtrack path pointer for all nodes in the 3d alignment cube for the 3 strings
	 * @param i - the current row position
	 * @param j - the current column position
	 * @param k - the current breadth position
	 * @param sbv - a buffer for constructing the first alignment string
	 * @param sbw - a buffer for constructing the second alignment string
	 * @param sbu - a buffer for constructing the third alignment string
	 * @param v - the first string to be aligned
	 * @param w - the second string to be aligned
	 * @param u - the third string to be aligned
	 */
	private static void ThreeDimensionalBacktrack(int[][][] backtrack, int i, int j, int k, StringBuffer sbv, StringBuffer sbw, StringBuffer sbu, String v, String w, String u) {
		if (i + j == 0) {
			sbv.insert(0, String.format("%" +(k) +"s", "").replace(' ', '-'));
			sbw.insert(0, String.format("%" +(k) +"s", "").replace(' ', '-'));
			sbu.insert(0, u.substring(0,k));
			return;
		}
		if (i + k == 0) {
			sbv.insert(0, String.format("%" +(j) +"s", "").replace(' ', '-'));
			sbw.insert(0, u.substring(0,j));
			sbu.insert(0, String.format("%" +(j) +"s", "").replace(' ', '-'));
			return;
		}
		if (j + k == 0) {
			sbv.insert(0, u.substring(0,i));
			sbw.insert(0, String.format("%" +(i) +"s", "").replace(' ', '-'));
			sbu.insert(0, String.format("%" +(i) +"s", "").replace(' ', '-'));
			return;
		}
		int pointer = backtrack[i][j][k];
		//System.out.println(pointer);
		if (pointer == ThreeDimensionalBacktrackPointers.INDEL_V.getId()) {
			sbv.append(v.charAt(i-1));
			sbw.append("-");
			sbu.append("-");
			ThreeDimensionalBacktrack(backtrack, i-1, j, k, sbv, sbw, sbu, v, w, u);
		} else if (pointer == ThreeDimensionalBacktrackPointers.INDEL_W.getId()) {
			sbv.append("-");
			sbw.append(w.charAt(j-1));
			sbu.append("-");
			ThreeDimensionalBacktrack(backtrack, i, j-1, k, sbv, sbw, sbu, v, w, u);
		} else if (pointer == ThreeDimensionalBacktrackPointers.INDEL_U.getId()) {
			sbv.append("-");
			sbw.append("-");
			sbu.append(u.charAt(k-1));
			ThreeDimensionalBacktrack(backtrack, i, j, k-1, sbv, sbw, sbu, v, w, u);
		} else if (pointer == ThreeDimensionalBacktrackPointers.INDEL_VW.getId()) {
			sbv.append(v.charAt(i-1));
			sbw.append(w.charAt(j-1));
			sbu.append("-");
			ThreeDimensionalBacktrack(backtrack, i-1, j-1, k, sbv, sbw, sbu, v, w, u);
		} else if (pointer == ThreeDimensionalBacktrackPointers.INDEL_VU.getId()) {
			sbv.append(v.charAt(i-1));
			sbw.append("-");
			sbu.append(u.charAt(k-1));
			ThreeDimensionalBacktrack(backtrack, i-1, j, k-1, sbv, sbw, sbu, v, w, u);
		} else if (pointer == ThreeDimensionalBacktrackPointers.INDEL_WU.getId()) {
			sbv.append("-");
			sbw.append(w.charAt(j-1));
			sbu.append(u.charAt(k-1));
			ThreeDimensionalBacktrack(backtrack, i, j-1, k-1, sbv, sbw, sbu, v, w, u);
		} else {
			sbv.append(v.charAt(i-1));
			sbw.append(w.charAt(j-1));
			sbu.append(u.charAt(k-1));
			ThreeDimensionalBacktrack(backtrack, i-1, j-1, k-1, sbv, sbw, sbu, v, w, u);
		}
	}
}

enum AlignmentType {
	GLOBAL(0), LOCAL(1), FITTING(2), OVERLAP(3);
	private final int id;
	AlignmentType(int id) {
		this.id = id;
	}
}

/**
 * a list of pointers determining the kind of movement/path that brings us from one node
 * to another (this is specific to matching 3 strings)
 */
enum ThreeDimensionalBacktrackPointers {
	INDEL_V(0), INDEL_W(1), INDEL_U(2), INDEL_VW(3), INDEL_VU(4), INDEL_WU(5), MATCH(6);
	private final int id;
	ThreeDimensionalBacktrackPointers(int id) {
		this.id = id;
	}
	
	public int getId() {
		return id;
	}
}

/**
 * an object used to represent the middle node and middle edge for the
 * divide and conquer algorithm for determining the alignment between
 * two strings
 */
class MiddleElement {
	int middleX;
	int middleY;
	int middleEdge;
	
	public MiddleElement(int middleX, int middleY, int middleEdge) {
		this.middleX = middleX;
		this.middleY = middleY;
		this.middleEdge = middleEdge;
	}
}

/**
 * an object used to represent the alignment between 2 strings 
 * and the alignment score for them
 */
class Alignment {
	int score = 0;
	String stringA, stringB;
	
	public Alignment(int score, String stringA, String stringB) {
		this.score = score;
		this.stringA = stringA;
		this.stringB = stringB;
	}
	
	@Override
	public String toString() {
		return score +"\r\n" +stringA +"\r\n" +stringB;
	}
	
	public int GetEditDistance() {
		int result = 0;
		for (int i = 0; i < stringA.length(); i++) {
			result += stringA.charAt(i) == stringB.charAt(i) ? 0 : 1;
		}
		return result;
	}
}

/**
 * represents one node in the path formed by a direct acylic graph
 * @variable predecessors - a list of other nodes whose outer edge points to this node
 * @variable weights - the weights of the edges in the path going to this node
 * @variable nodeIndex - an arbitrary id for the node (has no relation to the order of the node in the path) 
 * @variable backtrackPredecessor - a list of predecessor nodes that lead to the max score of this node
 */
class Node implements Comparable {
	List<Node> predecessors = new ArrayList<Node>();
	List<Integer> weights = new ArrayList<Integer>();
	int value = 0;
	int nodeIndex;
	List<Integer> backtrackPredecessor = new ArrayList<Integer>();
	boolean canBeAddedAsPredecessor = false;
	
	public Node(int nodeIndex) {
		this.nodeIndex = nodeIndex;
	}
	
	public void addPredecessor(Node n, int weight) {
		predecessors.add(n);
		weights.add(weight);
	}
	
	/**
	 * this method computes the highest score from source to this node while at the same time removing
	 * all other nodes that are not qualified to become predecessor node because they have no connection
	 * to the source node
	 * @param source - the first node in the path
	 */
	public void computeMaxValue(int source) {
		int prevValue, currValue;
		canBeAddedAsPredecessor = true;
		List<Integer> itemsForRemoval = new ArrayList<Integer>();
		for (int i = 0; i < predecessors.size(); i++) {
			if (predecessors.get(i).backtrackPredecessor.size() == 0 && 
					(predecessors.get(i).nodeIndex != source)) {
				itemsForRemoval.add(i);
			}
		}
		for (int i = itemsForRemoval.size() - 1; i >= 0 ; i--) {
			int removeIndex = itemsForRemoval.get(i);
			predecessors.remove(removeIndex);
			weights.remove(removeIndex);
		}
		for (int i = 0; i < predecessors.size(); i++) {
			prevValue = value;
			currValue = predecessors.get(i).value + weights.get(i);
			if (currValue > prevValue) {
				value = currValue;
				backtrackPredecessor = new ArrayList<Integer>();
				backtrackPredecessor.add(predecessors.get(i).nodeIndex);
			} else if (currValue == value) {
				backtrackPredecessor.add(predecessors.get(i).nodeIndex);
			}
		}
		List<Integer> temp = new ArrayList<Integer>();
		for (int i = 0; i < predecessors.size(); i++) {
			temp.add(predecessors.get(i).nodeIndex);
		}
	}
	
	@Override
	public String toString() {
		return "" +nodeIndex;
	}
	
	@Override
	public boolean equals(Object o){
		if(o instanceof Node){
		    Node toCompare = (Node) o;
		    return this.nodeIndex == toCompare.nodeIndex;
		}
		return false;
	}
	
	@Override
	public int compareTo(Object o) {  
		Node other = (Node)o;
		if (nodeIndex > other.nodeIndex)
			return 1;
		else if (nodeIndex < other.nodeIndex)
			return -1;
		return 0;
	}
}