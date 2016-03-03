import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.StringTokenizer;

public class ComputationalProteomics {
	
	public static void main(String[] args) {
		// common
		String path = "C:\\Users\\Timothy\\Desktop\\temp\\bioinformatics\\dataset_a_b.txt";
		List<String> lines = BioinformaticsCommon.readFile(path);
		
		// initialize shared variables
		List<SpectrumNode> spectrum = new ArrayList<SpectrumNode>();
		spectrum.add(new SpectrumNode(0));
		
		// test for constructing the graph of a spectrum
		/*List<SpectrumEdge> graph = ConstructSpectrumGraph(lines.get(0), spectrum);
		for (int i = 0; i < graph.size(); i++) {
			System.out.println(graph.get(i));
		}*/
		
		// test for generating ideal spectrum from peptide string
		//String temp = BioinformaticsCommon.joinList(IdealSpectrum("NQEL"), " ").trim();
		//System.out.println(temp);
		
		// test for decoding an ideal spectrum
		/*List<SpectrumEdge> graph = ConstructSpectrumGraph(lines.get(0), spectrum);
		String peptide = DecodingIdealSpectrum(graph, spectrum);
		System.out.println(peptide);*/
		
		// test for converting a peptide to a peptide vector
		//String temp = BioinformaticsCommon.joinList(ConvertToPeptideVector(lines.get(0)), " ");
		//BioinformaticsCommon.WriteOutputToFile(temp);
		
		// test for converting peptide vector into an amino acid string
		/*String temp = ConvertToPeptide(lines.get(0));
		System.out.println(temp);*/
		
		// test for peptide sequencing
		List<SpectrumEdge> edges = ConstructAllSimplePaths(lines.get(0), spectrum);
		String peptide = PeptideSequencing(spectrum, edges);
		System.out.println(peptide);
	}
	
	/**
	 * Finds all the amino acids that can be reconstructed from a given spectrum
	 * @param spectrum - a space delimeted list of integer masses representing peptide sequences
	 */
	private static List<SpectrumEdge> ConstructSpectrumGraph(String spectrumStr, List<SpectrumNode> spectrum) {
		StringTokenizer st = new StringTokenizer(spectrumStr);
		while (st.hasMoreTokens()) {
			spectrum.add(new SpectrumNode(Integer.parseInt(st.nextToken())));
		}
		List<SpectrumEdge> graph = new ArrayList<SpectrumEdge>();
		
		for (int i = 0; i < spectrum.size(); i++) {
			for (int j = i+1; j < spectrum.size(); j++) {
				int mass = spectrum.get(j).mass - spectrum.get(i).mass;
				if (BioinformaticsCommon.MASS_LIST_REV.containsKey(mass)) {
					graph.add(new SpectrumEdge(spectrum.get(i), spectrum.get(j), BioinformaticsCommon.MASS_LIST_REV.get(mass)));
				}
			}
		}
		
		return graph;
	}
	
	/**
	 * Generates the ideal spectrum from a string peptide
	 * @param peptide - an amino acid string
	 * @return
	 */
	private static List<Integer> IdealSpectrum(String peptide) {
		List<Integer> idealSpectrum = new ArrayList<Integer>();
		List<Integer> masses = new ArrayList<Integer>();
		
		masses.add(0);
		char[] acids = peptide.toCharArray();
		for (int i = 1; i <= acids.length; i++) {
			masses.add(i, masses.get(i-1) + BioinformaticsCommon.MASS_LIST.get(acids[i-1]+""));
		}
		
		idealSpectrum.add(0);
		for (int i = 0; i < masses.size(); i++) {
			for (int j = i+1; j < masses.size(); j++) {
				idealSpectrum.add(masses.get(j) - masses.get(i));
			}
		}
		Collections.sort(idealSpectrum);
		return idealSpectrum;
	}
	
	private static void DFS(List<SpectrumEdge> graph, int currNode, int targetNode, String currStr, List<String> strList) {
		boolean backtrackMarker = false;
		for (int i = 0; i < graph.size(); i++) {
			if (graph.get(i).nodeA.mass == currNode) {
				if (backtrackMarker) {
					currStr = currStr.substring(0, currStr.length()-1);
				} else {
					backtrackMarker = true;
				}
				currStr += graph.get(i).aminoAcid;
				if (graph.get(i).nodeB.mass == targetNode) {
					strList.add(currStr);
				} else {
					DFS(graph, graph.get(i).nodeB.mass, targetNode, currStr, strList);
				}
			}
			
		}
	}
	
	/**
	 * Generates the peptide string that explains the spectrum
	 * @param graph - the list of edges generated from the spectrum
	 * @param spectrum - a list of integer masses
	 * @return
	 */
	private static String DecodingIdealSpectrum(List<SpectrumEdge> graph, List<SpectrumNode> spectrum) {
		List<String> strList = new ArrayList<String>();
		DFS(graph, 0, graph.get(graph.size()-1).nodeB.mass, "", strList);
		for (int i = 0; i < strList.size(); i++) {
			List<Integer> idealSpectrum = IdealSpectrum(strList.get(i));
			boolean isSubset = true;
			for (int j = 0; j < spectrum.size(); j++) {
				if (!idealSpectrum.contains(spectrum.get(j).mass)) {
					isSubset = false;
					break;
				}
			}
			if (isSubset) {
				return strList.get(i);
			}
		}
		return "no matching ideal spectrum";
	}
	
	/**
	 * Converts a peptide into a peptide vector
	 * @param peptide - an amino acid string
	 * @return a list of integers containing 1 at each of the prefix coordinates (the array index pertaining to the mass of the peptide) and 0 otherwise
	 */
	private static List<Integer> ConvertToPeptideVector(String peptide) {
		List<Integer> peptideVector = new ArrayList<Integer>();
		for (int i = 0; i < peptide.length(); i++) {
			String aminoAcid = peptide.substring(i, i+1);
			int mass = BioinformaticsCommon.MASS_LIST.get(aminoAcid);
			for (int j = 0; j < mass; j++) {
				if (j == mass-1) {
					peptideVector.add(1);
				} else {
					peptideVector.add(0);
				}
			}
		}
		return peptideVector;
	}
	
	/**
	 * Converts a peptide vector into an amino acid string
	 * @param peptideVector - a list of integers containing 1 at each of the prefix coordinates (the array index pertaining to the mass of the peptide) and 0 otherwise
	 * @return the amino acid string that represents this peptide vector
	 */
	private static String ConvertToPeptide(String peptideVector) {
		String peptide = "";
		StringTokenizer st = new StringTokenizer(peptideVector);
		int index = 0;
		while (st.hasMoreTokens()) {
			index++;
			String token = st.nextToken();
			if ("1".equals(token)) {
				peptide += BioinformaticsCommon.MASS_LIST_REV.get(index);
				index = 0;
			}
		}
		return peptide;
	}
	
	/**
	 * Connect all pairs of nodes with index difference that is in the list of integer protein masses
	 * @param spectrumStr - a space delimeted spectral vector
	 * @param spectrum - a storage list for each score in spectrumStr; already contains 0 as the first element
	 * @return
	 */
	private static List<SpectrumEdge> ConstructAllSimplePaths(String spectrumStr, List<SpectrumNode> spectrum) {
		StringTokenizer st = new StringTokenizer(spectrumStr);
		int index = 0;
		while (st.hasMoreTokens()) {
			spectrum.add(new SpectrumNode(Integer.parseInt(st.nextToken()), ++index));
		}
		List<SpectrumEdge> graph = new ArrayList<SpectrumEdge>();
		
		for (int i = 0; i < spectrum.size(); i++) {
			for (int j = i+1; j < spectrum.size(); j++) {
				int nodeDiff = j - i;
				if (BioinformaticsCommon.MASS_LIST_REV.containsKey(nodeDiff)) {
					graph.add(new SpectrumEdge(spectrum.get(i), spectrum.get(j), BioinformaticsCommon.MASS_LIST_REV.get(nodeDiff)));
				}
			}
		}
		return graph;
	}
	
	/**
	 * Using the Bellman-Ford algorithm to solve for the longest path (by negating the mass), get an amino acid string that maximizes 
	 * the score against an input spectral vector
	 * @param vertices - space delimited spectral vector
	 * @param edges - all the edges formed by connecting a pair of nodes with an index difference that is equivalent to any of the protein masses
	 * @return an amino acid string that maximizes the score against the spectral vector
	 */
	private static String PeptideSequencing(List<SpectrumNode> vertices, List<SpectrumEdge> edges) {
		int size = vertices.size();
		int[] distance = new int[size];
		int[] predecessor = new int[size];
		final int INFINITY = 99999; //do not use Integer.MAX_VALUE
		
		for (int v = 0; v < size; v++) {
			distance[v] = INFINITY;
			predecessor[v] = -1;
		}
		
		distance[0] = 0;
		
		for (int i = 0; i < size; i++) {
			for (int j = 0; j < edges.size(); j++) {
				SpectrumEdge edge = edges.get(j);
				int u = edge.nodeA.index;
				int v = edge.nodeB.index;
				int w = -1 * edge.nodeB.mass;
				if (/*u == i && */distance[u] + w < distance[v]) {
					if (distance[u] == INFINITY)
						distance[v] = INFINITY;
					else
						distance[v] = distance[u] + w;
					predecessor[v] = u;
				}
			}
		}
		
		
		int lastIndexChecked = size-1;
		String peptide = "";
		int score = 0;
		while (lastIndexChecked != 0) {
			int mass = lastIndexChecked - predecessor[lastIndexChecked];
			peptide = BioinformaticsCommon.MASS_LIST_REV.get(mass) + peptide;
			//System.out.println(predecessor[lastIndexChecked] +" ~ " +lastIndexChecked +":" +distance[lastIndexChecked]);
			score += distance[lastIndexChecked];
			lastIndexChecked = predecessor[lastIndexChecked];
		}
		return peptide;
	}
	
}

class SpectrumNode {
	int mass;
	int index;
	
	public SpectrumNode(int mass) {
		this.mass = mass;
	}
	
	public SpectrumNode(int mass, int index) {
		this.mass = mass;
		this.index = index;
	}
}

class SpectrumEdge {
	SpectrumNode nodeA;
	SpectrumNode nodeB;
	String aminoAcid;
	
	public SpectrumEdge(SpectrumNode nodeA, SpectrumNode nodeB, String aminoAcid) {
		this.nodeA = nodeA;
		this.nodeB = nodeB;
		this.aminoAcid = aminoAcid;
	}
	
	@Override
	public String toString() {
		return this.nodeA.mass +BioinformaticsCommon.NODE_SEPARATOR +this.nodeB.mass +":" +this.aminoAcid;
	}
}