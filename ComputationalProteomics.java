import java.util.ArrayList;
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
		List<SpectrumEdge> graph = ConstructSpectrumGraph(lines.get(0), spectrum);
		for (int i = 0; i < graph.size(); i++) {
			System.out.println(graph.get(i));
		}	
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
}

class SpectrumNode {
	int mass;
	
	public SpectrumNode(int mass) {
		this.mass = mass;
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