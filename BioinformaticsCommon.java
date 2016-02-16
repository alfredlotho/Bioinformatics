import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.StringTokenizer;

public class BioinformaticsCommon {
	
	public static final String NODE_SEPARATOR = "->";

	/**
	 * reads a file into memory
	 * @return a list of strings in the file (line by line)
	 */
	public static List<String> readFile() {
		String path = "C:\\Users\\Alfred\\Downloads\\dataset_a_b.txt";
		return readFile(path);
	}
	
	public static List<String> readFile(String path) {
		String line = "";
		List<String> records = new ArrayList<String>();
		
		BufferedReader br;
		
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(path), "UTF-8"));
			int i = 0;
			line = br.readLine();
			while (line != null) {
				records.add(line);
				line = br.readLine();				
			}
			br.close();
		} catch (Exception e) {
			
		} 
		
		return records;
	}
	
	/**
	 * @param the nucleotide string you want to perform a reverse complement on
	 * @return the reverse complement
	 * Sample Input:
	 * CGGCT
	 * Sample Output:
	 * AGCCG
	 */
	public static String ReverseComplement(String text) {
		String reverseComplement = "";
		String forward = "CGAT";
		String reverse = "GCTA";
		for (int i = text.length()-1; i >= 0; i--) {
			int matchingIndex = forward.indexOf(text.charAt(i));
			reverseComplement += reverse.charAt(matchingIndex);
		}
		return reverseComplement;
	}
	
	/**
	 * concatenates all the nucleotides/peptides in a list
	 */
	public static <T> String joinPath(List<T> stringList) {
		String retString = stringList.get(0).toString();
		for (int i = 1; i < stringList.size(); i++) {
			String currString = stringList.get(i).toString();
			retString += currString.charAt(currString.length()-1);
		}
		return retString;
	}
	
	/**
	 * concatenates a pair of nucleotides/peptides in a list
	 * @param stringList - contains a list of strings, each having the following format: 
	 * <nucleotide1>|<nucleotide2> (e.g. A|G)
	 * @return - 2 nucleotide strings formed from the first and second tokens of the strings, delimeted
	 * by the pipe character (|), in the list
	 * Sample Input:
	 * {"A|G", "C|G", "T|T"}
	 * Sample Output:
	 * "ACT"
	 * "GGT"
	 */
	public static <T> List<String> joinPathPair(List<T> stringList) {
		StringTokenizer st = new StringTokenizer(stringList.get(0).toString(), "|");
		String pairOne = st.nextToken();
		String pairTwo = st.nextToken();
		for (int i = 1; i < stringList.size(); i++) {
			st = new StringTokenizer(stringList.get(i).toString(), "|");
			String tokenOne = st.nextToken();
			String tokenTwo = st.nextToken();
			pairOne += tokenOne.charAt(tokenOne.length()-1);
			pairTwo += tokenTwo.charAt(tokenTwo.length()-1);
		}
		
		List<String> retStringList = new ArrayList<String>();
		retStringList.add(pairOne);
		retStringList.add(pairTwo);

		return retStringList;
	}
	
	/**
	 * concatenates all strings in the list and separates them with the separator character
	 */
	public static <T> String joinList(List<T> stringList, String separator) {
		String retString = "";
		for (int i = 0; i < stringList.size(); i++) {
			if (i != 0)
				retString += separator;
			retString += stringList.get(i);
		}
		return retString;
	}
	
	/**
	 * concatenates all strings in the list in reverse order and separates them with the 
	 * separator character
	 */
	public static <T> String joinListReverse(List<T> stringList, String separator) {
		String retString = "";
		for (int i = stringList.size() - 1; i >= 0 ; i--) {
			if (i != stringList.size() - 1)
				retString += separator;
			retString += stringList.get(i);
		}
		return retString;
	}
	
	/**
	 * concatenates all strings in the list and encloses them in a pair of parentheses
	 */
	public static <T> String joinListWithEnclosure(List<T> stringList) {
		String retString = "";
		for (int i = 0; i < stringList.size(); i++) {
			retString += "(" +stringList.get(i) +")";
		}
		return retString;
	}
	
	/**
	 * joins all the peptides in the list, converts them into their mass equivalents and
	 * separates them with a dash character
	 * Sample Input:
	 * {"G", "A", "S"}
	 * Sample Output:
	 * 57-71-87
	 */
	public static String joinPeptidesAndConvertToMass(List<String> stringList) {
		String retString = "";
		for (int i = 0; i < stringList.size(); i++) {
			if (i != 0)
				retString += " ";
			//retString += stringList.get(i);
			String peptide = stringList.get(i);
			for (int j = 0; j < peptide.length(); j++) {
				String aminoAcid = peptide.charAt(j) + "";
				if (j != 0) 
					retString += "-";
				retString += MASS_LIST.get(aminoAcid);
			}
		}
		return retString;
	}
	
	/**
	 * converts the peptide to its mass equivalent
	 */
	public static String convertToMass(String peptide) {
		String retString = "";
		for (int j = 0; j < peptide.length(); j++) {
			String aminoAcid = peptide.charAt(j) + "";
			if (j != 0) 
				retString += "-";
			retString += MASS_LIST.get(aminoAcid);
		}
		return retString;
	}
	
	/**
	 * concatenates all strings in the set and separates them with the separator character
	 */
	public static <T> String joinSet(Set<T> stringSet, String separator) {
		String retString = "";
		Iterator iter = stringSet.iterator();
		while (iter.hasNext()) {
			retString += iter.next();
			if (iter.hasNext())
				retString += separator;
		}
		
		return retString;
	}
	
	/**
	 * writes a string to a text file
	 */
	public static void WriteOutputToFile(String temp) {
		String path = "C:\\Users\\Timothy\\Downloads\\output.txt";
		try (Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(path), "utf-8"))) {
		   writer.write(temp);
		} catch (Exception e) {
			
		}
	}
	
	public static int RandomWithRange(int min, int max)
	{
	   int range = (max - min) + 1;     
	   return (int)(Math.random() * range) + min;
	}
	
	public static final Map<String, String> CODON_LIST;
	static
	{
		CODON_LIST = new HashMap<String, String>();
		CODON_LIST.put("AAA","K");
		CODON_LIST.put("AAC","N");
		CODON_LIST.put("AAG","K");
		CODON_LIST.put("AAU","N");
		CODON_LIST.put("ACA","T");
		CODON_LIST.put("ACC","T");
		CODON_LIST.put("ACG","T");
		CODON_LIST.put("ACU","T");
		CODON_LIST.put("AGA","R");
		CODON_LIST.put("AGC","S");
		CODON_LIST.put("AGG","R");
		CODON_LIST.put("AGU","S");
		CODON_LIST.put("AUA","I");
		CODON_LIST.put("AUC","I");
		CODON_LIST.put("AUG","M");
		CODON_LIST.put("AUU","I");
		CODON_LIST.put("CAA","Q");
		CODON_LIST.put("CAC","H");
		CODON_LIST.put("CAG","Q");
		CODON_LIST.put("CAU","H");
		CODON_LIST.put("CCA","P");
		CODON_LIST.put("CCC","P");
		CODON_LIST.put("CCG","P");
		CODON_LIST.put("CCU","P");
		CODON_LIST.put("CGA","R");
		CODON_LIST.put("CGC","R");
		CODON_LIST.put("CGG","R");
		CODON_LIST.put("CGU","R");
		CODON_LIST.put("CUA","L");
		CODON_LIST.put("CUC","L");
		CODON_LIST.put("CUG","L");
		CODON_LIST.put("CUU","L");
		CODON_LIST.put("GAA","E");
		CODON_LIST.put("GAC","D");
		CODON_LIST.put("GAG","E");
		CODON_LIST.put("GAU","D");
		CODON_LIST.put("GCA","A");
		CODON_LIST.put("GCC","A");
		CODON_LIST.put("GCG","A");
		CODON_LIST.put("GCU","A");
		CODON_LIST.put("GGA","G");
		CODON_LIST.put("GGC","G");
		CODON_LIST.put("GGG","G");
		CODON_LIST.put("GGU","G");
		CODON_LIST.put("GUA","V");
		CODON_LIST.put("GUC","V");
		CODON_LIST.put("GUG","V");
		CODON_LIST.put("GUU","V");
		CODON_LIST.put("UAA","" );
		CODON_LIST.put("UAC","Y");
		CODON_LIST.put("UAG","" );
		CODON_LIST.put("UAU","Y");
		CODON_LIST.put("UCA","S");
		CODON_LIST.put("UCC","S");
		CODON_LIST.put("UCG","S");
		CODON_LIST.put("UCU","S");
		CODON_LIST.put("UGA","" );
		CODON_LIST.put("UGC","C");
		CODON_LIST.put("UGG","W");
		CODON_LIST.put("UGU","C");
		CODON_LIST.put("UUA","L");
		CODON_LIST.put("UUC","F");
		CODON_LIST.put("UUG","L");
		CODON_LIST.put("UUU","F");
	}
	
	public static final Map<String, Integer> MASS_LIST;
	static
	{
		MASS_LIST = new HashMap<String, Integer>();
		MASS_LIST.put("G",57);
		MASS_LIST.put("A",71);
		MASS_LIST.put("S",87);
		MASS_LIST.put("P",97);
		MASS_LIST.put("V",99);
		MASS_LIST.put("T",101);
		MASS_LIST.put("C",103);
		MASS_LIST.put("I",113);
		MASS_LIST.put("L",113);
		MASS_LIST.put("N",114);
		MASS_LIST.put("D",115);
		MASS_LIST.put("K",128);
		MASS_LIST.put("Q",128);
		MASS_LIST.put("E",129);
		MASS_LIST.put("M",131);
		MASS_LIST.put("H",137);
		MASS_LIST.put("F",147);
		MASS_LIST.put("R",156);
		MASS_LIST.put("Y",163);
		MASS_LIST.put("W",186);
	}
	
	public static final Map<String, Integer> MASS_LIST_TRIM;
	static
	{
		MASS_LIST_TRIM = new LinkedHashMap<String, Integer>();
		MASS_LIST_TRIM.put("G",57);
		MASS_LIST_TRIM.put("A",71);
		MASS_LIST_TRIM.put("S",87);
		MASS_LIST_TRIM.put("P",97);
		MASS_LIST_TRIM.put("V",99);
		MASS_LIST_TRIM.put("T",101);
		MASS_LIST_TRIM.put("C",103);
		MASS_LIST_TRIM.put("I",113);
		MASS_LIST_TRIM.put("N",114);
		MASS_LIST_TRIM.put("D",115);
		MASS_LIST_TRIM.put("K",128);
		MASS_LIST_TRIM.put("E",129);
		MASS_LIST_TRIM.put("M",131);
		MASS_LIST_TRIM.put("H",137);
		MASS_LIST_TRIM.put("F",147);
		MASS_LIST_TRIM.put("R",156);
		MASS_LIST_TRIM.put("Y",163);
		MASS_LIST_TRIM.put("W",186);
	}
	
	private static String scoringMatrixPath = "C:\\Users\\Alfred\\Downloads\\";
	private static String scoringMatrixType = "BLOSUM62";
	public static void SetScoringMatrix(String id) {
		scoringMatrixType = id;
	}
	
	public static List<String> AMINO_ACID_LIST_FOR_SCORING = new ArrayList<String>();
	
	/**
	 * loads the scoring matrix (e.g. BLOSUM62, PAM250) form a text file to a 2d array
	 * @return
	 */
	public static int[][] GetScoringMatrix() {
		String line = "";
		int[][] scoringMatrix = null;
		String path = scoringMatrixPath +scoringMatrixType + ".txt";
		BufferedReader br;
		
		try {
			br = new BufferedReader(new InputStreamReader(new FileInputStream(path), "UTF-8"));
			line = br.readLine();
			StringTokenizer st = new StringTokenizer(line);
			while(st.hasMoreTokens()) {
				AMINO_ACID_LIST_FOR_SCORING.add(st.nextToken());
			}
			int matrixSize = AMINO_ACID_LIST_FOR_SCORING.size();
			scoringMatrix = new int[matrixSize][matrixSize];
			line = br.readLine();
			int i = 0;
			while (line != null) {
				st = new StringTokenizer(line.substring(2, line.length()));
				int j = 0;
				while (st.hasMoreTokens()) {
					scoringMatrix[i][j++] = Integer.parseInt(st.nextToken());
				}
				i++;
				line = br.readLine();
			}
			br.close();
		} catch (Exception e) {
			
		} 
		
		return scoringMatrix;
	}
	
	/**
	 * @param s - a dna string or its numeric node index 
	 * 			- we are sure that this is either a positive number of a nonempty string so no need for those checks
	 * @return true if the input string is an integer
	 */
	public static boolean isInteger(String s) {
	    for(int i = 0; i < s.length(); i++) {
	        if(Character.digit(s.charAt(i), 10) < 0) {
	        	return false;
	        }
	    }
	    return true;
	}

}