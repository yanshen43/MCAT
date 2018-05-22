package helper;
/*NOTE: Currently the SelectAndRemoveForTesting only uses the last 2/3 for training and the first 1/3 for testing.
 * (Instead of randomly selection)
 */


import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Scanner;
import java.util.Set;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;
import org.apache.commons.math.stat.ranking.NaNStrategy;
import org.apache.commons.math.stat.ranking.NaturalRanking;
import org.apache.commons.math.stat.ranking.TiesStrategy;



public class SeqTools {
	
	private static double SCORING_CUTOFF_PERCENTAGE = 0.7;
	public static String convertToDigitRepresentation(String seq) {
		seq = seq.replace('A', '0');
		seq = seq.replace('C', '1');
		seq = seq.replace('G', '2');
		seq = seq.replace('T', '3');
		return seq;
	}

	public static String convertToDNARepresentation(String seq) {
		seq = seq.replace('0', 'A');
		seq = seq.replace('1', 'C');
		seq = seq.replace('2', 'G');
		seq = seq.replace('3', 'T');
		return seq;
	}
	public static void normalize(Hashtable<String, double[]> hash) {
		double sum0 = 0, sum1 = 0;
		Enumeration<String> e = hash.keys();
		while (e.hasMoreElements()) {
			String nextKey = e.nextElement();
			sum0 += (double)hash.get(nextKey)[0];
			sum1 += (double)hash.get(nextKey)[1];
		}
		e = hash.keys();
		while (e.hasMoreElements()) {
			String nextKey = e.nextElement();
			double[] freqs = hash.get(nextKey);
			freqs[0] = freqs[0] / sum0;
			freqs[1] = freqs[1] / sum1;
			freqs[2] = freqs[1]-freqs[0];
			hash.put(nextKey, freqs);
		}
	}

	/*	public static Hashtable<String, Double> getKMerComposition(ArrayList<String> seqs, int MOTIF_WIDTH, boolean limitKMerCountPerSequenceTo2, boolean ignoreSimpleRepeats) {
		Hashtable<String, Double> kMerFreqs = new Hashtable<String, Double>();
		for (String s : seqs) {
			int len = s.length();
			Hashtable<String, Integer> kMerCountsInThisSeq = new Hashtable<String, Integer>();
			for (int i = 0; i < s.length()-MOTIF_WIDTH+1; i++) {
				String substr = s.substring(i, i+MOTIF_WIDTH);
				if (substr.contains("N")) {
					continue;
				}
				if (ignoreSimpleRepeats) {
					//Ignore simple repeats
					//If >=75% positions in the substring is one nucleotide, this is considered a simple repeat
					int[] baseCounter = new int[4];
					for (int j = 0; j < substr.length(); j++) {
						baseCounter[substr.charAt(j)-48]++;
					}
					boolean repeat = false;
					for (int counter:baseCounter) {
						if (counter >= Math.ceil(MOTIF_WIDTH * 0.75)) {
							repeat = true;
						}
					}
					int counterZero = 0;
					for (int counter:baseCounter) {
						if (counter==0) {
							counterZero++;
						}
					}
					if (counterZero>=2) {
						for (int counter:baseCounter) {
							if (counter >= Math.floor(MOTIF_WIDTH * 0.5)) {
								repeat = true;
							}
						}

					}
					//Remove dinucleotide repeats
					Hashtable<String, Integer> dinucleotideCounts = new Hashtable<String, Integer>();
					for (int j = 0; j < substr.length()-1; j++) {
						String dinucleotide = substr.substring(j, j+2);
						if (dinucleotideCounts.containsKey(dinucleotide)) {
							dinucleotideCounts.put(dinucleotide, dinucleotideCounts.get(dinucleotide)+1);
						} else {
							dinucleotideCounts.put(dinucleotide, 1);
						}
					}

					for (String kMer: dinucleotideCounts.keySet()) {
						int count = dinucleotideCounts.get(kMer);
						if (count >= (int) MOTIF_WIDTH /2) {
							repeat = true;
						}
					}
					if (repeat == true) {
						continue;
					}
				}
				if (kMerCountsInThisSeq.containsKey(substr)) {
					kMerCountsInThisSeq.put(substr, kMerCountsInThisSeq.get(substr)+1);
				} else {
					kMerCountsInThisSeq.put(substr, 1);
				}
			}
			Enumeration<String> e = kMerCountsInThisSeq.keys();
			while (e.hasMoreElements()) {
				String kMer = e.nextElement();
				double freq;
				if ((limitKMerCountPerSequenceTo2) && (kMerCountsInThisSeq.get(kMer)>=2)) {
					freq = 2.0/(double)len;
				}
				else {
					freq = (double)kMerCountsInThisSeq.get(kMer)/(double)len;
				}
				if (kMerFreqs.containsKey(kMer)) {
					kMerFreqs.put(kMer, kMerFreqs.get(kMer) + freq);
				}
				else {
					kMerFreqs.put(kMer, freq);
				}
			}
		}
		return kMerFreqs;
	}*/

	public static boolean isSimpleRepeat(String kMer) {
		//Ignore simple repeats
		//If >75% positions in the substring is one nucleotide, this is considered a simple repeat
		int[] baseCounter = new int[4];
		for (int j = 0; j < kMer.length(); j++) {
			baseCounter[kMer.charAt(j)-48]++;
		}
		boolean isRepeat = false;
		for (int counter:baseCounter) {
			if ((kMer.length()>=8) && (counter >= Math.floor(kMer.length() * 0.75))) {
				isRepeat = true;
			}
		}
		int counterZero = 0; //count how many bases (of the four) exist in this kmer
		for (int counter:baseCounter) {
			if (counter==0) {
				counterZero++;
			}
		}
		if ((counterZero >= 2) && (kMer.length()>=8)) {
			isRepeat = true;
		}
		
		
		//If one base occurs continuously for >=50% positions then this is considered to be a simple repeat 
		if (!isRepeat && kMer.length()>=8) {				
			int continuousSingleBase = 1;
			char previousBase = kMer.charAt(0);
			for (int j = 1; j < kMer.length(); j++) {
				if (kMer.charAt(j) == previousBase) {
					continuousSingleBase++;
					if (continuousSingleBase >= Math.floor(kMer.length()/2)) {
						isRepeat = true;
						break;
					}
				}
				else {
					continuousSingleBase = 1;
					previousBase = kMer.charAt(j);
				}
			}


		}
		return isRepeat;
	}
	public static boolean isDinucleotideRepeat(String kMer) {
		//Remove dinucleotide repeats
		boolean isRepeat = false;
		int[] baseCounter = new int[4];
		for (int j = 0; j < kMer.length(); j++) {
			baseCounter[kMer.charAt(j)-48]++;
		}
		int counterZero = 0; //count how many bases (of the four) exist in this kmer
		for (int counter:baseCounter) {
			if (counter==0) {
				counterZero++;
			}
		}
		
		
		for (int phase = 0; phase <= 1; phase++) {
			
			int counterDinucleotideRepeats  = 0;
			String previousDinucleotide = null;
			for (int j = phase; j < kMer.length()-1; j=j+2) {
				String dinucleotide = kMer.substring(j, j+2);
				
				if (previousDinucleotide == null) {
					previousDinucleotide = dinucleotide;
					counterDinucleotideRepeats++;
				}
				else {
					if (dinucleotide.equals(previousDinucleotide)) {
						counterDinucleotideRepeats++;
					} else {
						previousDinucleotide = dinucleotide;
						counterDinucleotideRepeats = 1;
					}
				}
			}
			
			
			Hashtable<String, Integer> dinucleotideCounts = new Hashtable<String, Integer>();
			for (int j = phase; j < kMer.length()-1; j=j+2) {
				String dinucleotide = kMer.substring(j, j+2);
				if (dinucleotideCounts.containsKey(dinucleotide)) {
					dinucleotideCounts.put(dinucleotide, dinucleotideCounts.get(dinucleotide)+1);
				} else {
					dinucleotideCounts.put(dinucleotide, 1);
				}
			}
			
			int maxIndividualDinucleotideCount = 0;
			for (String kmer1: dinucleotideCounts.keySet()) {
				if (dinucleotideCounts.get(kmer1) > maxIndividualDinucleotideCount) {
					maxIndividualDinucleotideCount = dinucleotideCounts.get(kmer1);
				}
			}
			/*System.out.println("counterZero="+counterZero);
			System.out.println("max individual dinucleotideCount = "+maxIndividualDinucleotideCount);*/
			
			if ((counterDinucleotideRepeats > Math.ceil(kMer.length() /3)) || ((counterZero>=2) && maxIndividualDinucleotideCount >= Math.floor(kMer.length() /3))) {
				isRepeat = true;
				break;
			}
		}
		return isRepeat;
	}



//Positive: 1; negative: 0;
//the double[] in kMerFreqs: 0 - in negative sequences; 1 - in positive sequences; 2 - frequency difference (pos-neg), generated after calling normalize
public static void getKMerComposition(ArrayList<String> seqs, int MOTIF_WIDTH, boolean limitKMerCountPerSequence, int upperLimitOfKMerCountPerSeq, Hashtable<String, double[]> kMerFreqs, int isPositive, boolean ignoreSimpleRepeats) {
	for (String s : seqs) {
		int len = s.length();
		Hashtable<String, Integer> kMerCountsInThisSeq = new Hashtable<String, Integer>();
		for (int i = 0; i < s.length()-MOTIF_WIDTH+1; i++) {
			String substr = s.substring(i, i+MOTIF_WIDTH);
			if (substr.contains("N")) {
				continue;
			}
			if (ignoreSimpleRepeats) {
				if (isSimpleRepeat(substr) || isDinucleotideRepeat(substr)) {
					continue;
				}
			}

			if (kMerCountsInThisSeq.containsKey(substr)) {
				kMerCountsInThisSeq.put(substr, kMerCountsInThisSeq.get(substr)+1);
			} else {
				kMerCountsInThisSeq.put(substr, 1);
			}
		}
		Enumeration<String> e = kMerCountsInThisSeq.keys();
		while (e.hasMoreElements()) {
			String kMer = e.nextElement();
			double freq;
			if ((limitKMerCountPerSequence) && (kMerCountsInThisSeq.get(kMer)>=upperLimitOfKMerCountPerSeq)) {
				freq = (double)upperLimitOfKMerCountPerSeq / (double)len;
			}
			else {
				freq = (double)kMerCountsInThisSeq.get(kMer) / (double)len;
			}
			if (kMerFreqs.containsKey(kMer)) {
				double[] counts = kMerFreqs.get(kMer);
				counts[isPositive] += freq;
				kMerFreqs.put(kMer, counts);
			}
			else {
				double[] counts = new double[3];
				counts[isPositive] = freq;
				kMerFreqs.put(kMer, counts);
			}
		}
	}
}

public static Pwm calculateBackground(ArrayList<String> seqs) {
	int[] counter = new int[4];
	for (String s : seqs) {
		char[] letters = s.toCharArray();
		for (char c : letters) {
			if (c=='N') {
				continue;
			}
			counter[c-48]++;
		}
	}
	int sum = 0;
	for (int i = 0; i < 4; i++) {
		sum += counter[i];
	}
	double[][] matrix = new double[4][1];
	for (int i=0; i<4; i++) {
		matrix[i][0] = (double)counter[i] / (double)sum;
	}
	Pwm background = new Pwm(matrix);
	return background;
}

public static double[][] getLogScoringMatrix(Pwm theta, Pwm background) {
	double[][] logScoringMatrix = theta.getMatrix();
	double[][] backgroundMatrix = background.getMatrix();
	int motifWidth = logScoringMatrix[0].length;
	for (int i = 0; i < 4; i++ ) {
		for (int j = 0; j < motifWidth; j++) {
			if (logScoringMatrix[i][j]==0) {
				logScoringMatrix[i][j] = 0.001;
			}
			if (backgroundMatrix[i][0]==0) {
				backgroundMatrix[i][0] = 0.001;
			}
			logScoringMatrix[i][j] = Math.log(logScoringMatrix[i][j]) - Math.log(backgroundMatrix[i][0]);
		}
	}

	return logScoringMatrix;
}
public static double[] getMaxScores(ArrayList<String> seqs, Pwm theta, Pwm background) {
	int size = seqs.size();
	double[] scores = new double[size];
	double[][] logScoringMatrix = getLogScoringMatrix(theta, background);
	int i = 0;
	for (String s : seqs) {
		MotifInstance r = scanMotifMax(s, logScoringMatrix);
		if (r != null) {
			scores[i] = r.getScore();
		} else {
			scores[i] = 0;
		}
		i++;
	}
	return scores;
}


public static double rankSumTest(ArrayList<String> posSeqs, ArrayList<String> negSeqs, Pwm theta, Pwm bg) throws MathException {
	double[] maxScoresPos = getMaxScores(posSeqs, theta, bg);
	double[] maxScoresNeg = getMaxScores(negSeqs, theta, bg);

	int n1 = maxScoresPos.length;
	int n2 = maxScoresNeg.length;

	double[] maxScoresAll = new double[n1 + n2];

	for (int i = 0; i < n1; i++) {
		maxScoresAll[i] = maxScoresPos[i];
	}
	for (int i = 0; i < n2; i++) {
		maxScoresAll[i + n1] = maxScoresNeg[i];
	}

	NaturalRanking ranking = new NaturalRanking(NaNStrategy.REMOVED, TiesStrategy.AVERAGE);
	double[] ranks = ranking.rank(maxScoresAll);

	double tStatistic = 0.0;
	for (int i = 0; i < n1; i++) {
		tStatistic += ranks[i];
	}

	double tMean = (double)n1 * (n1 + n2 + 1) / 2.0;
	double tStd = Math.sqrt((double)n1 * n2 / 12.0 * (n1 + n2 + 1));
	double zStatistic = (tStatistic - tMean) / tStd;

	NormalDistribution standardNormal = new NormalDistributionImpl(0, 1);
	double pvalue = 1.0 - standardNormal.cumulativeProbability(zStatistic);

	return pvalue;
}

public static double rankSumTest(ArrayList<String> posSeqs, ArrayList<String> negSeqs, Pwm theta, Pwm bgPos, Pwm bgNeg) throws MathException {
	double[] maxScoresPos = getMaxScores(posSeqs, theta, bgPos);
	double[] maxScoresNeg = getMaxScores(negSeqs, theta, bgNeg);

	int n1 = maxScoresPos.length;	int n2 = maxScoresNeg.length;

	double[] maxScoresAll = new double[n1 + n2];

	for (int i = 0; i < n1; i++) {
		maxScoresAll[i] = maxScoresPos[i];
	}
	for (int i = 0; i < n2; i++) {
		maxScoresAll[i + n1] = maxScoresNeg[i];
	}

	NaturalRanking ranking = new NaturalRanking(NaNStrategy.REMOVED, TiesStrategy.AVERAGE);
	double[] ranks = ranking.rank(maxScoresAll);

	double tStatistic = 0.0;
	for (int i = 0; i < n1; i++) {
		tStatistic += ranks[i];
	}

	double tMean = (double)n1 * (n1 + n2 + 1) / 2.0;
	double tStd = Math.sqrt((double)n1 * n2 / 12.0 * (n1 + n2 + 1));
	double zStatistic = (tStatistic - tMean) / tStd;

	NormalDistribution standardNormal = new NormalDistributionImpl(0, 1);
	double pvalue = 1.0 - standardNormal.cumulativeProbability(zStatistic);

	return pvalue;
}

public static FastaSequences loadSeqs(String filename, boolean checkBothStrands) throws FileNotFoundException {
	ArrayList<String> seqs = new ArrayList<String>();
	ArrayList<String> ids = new ArrayList<String>();
	Scanner inputStream = null;
	String line = "";

	inputStream=new Scanner(new FileInputStream(filename));  

	//line = inputStream.nextLine().trim().toUpperCase();
	while (inputStream.hasNextLine()) {
		line = inputStream.nextLine().trim();			
		if (line.startsWith(">")) {
			String id = line.substring(1);
			ids.add(id);
			if (checkBothStrands) {
				ids.add(id+"|revcom");
			}
			continue;
		}
		if (line.isEmpty()) {
			continue;
		}

		line = line.toUpperCase();
		StringBuffer sb = new StringBuffer();
		
		do {
			if (line.contains("A") || line.contains("T") || line.contains("C") || line.contains("G")) { //Convert the nucleotides to numbers 0-3
				char[] arr = line.toCharArray();
				for (char i : arr) {
					switch (i) {
					case 'A': sb.append("0");break;
					case 'C': sb.append("1");break;
					case 'G': sb.append("2");break;
					case 'T': sb.append("3");break;
					default: sb.append("N");break;
					}
				}
			}
			if (inputStream.hasNextLine()) {
				line = inputStream.nextLine().trim();
				if (line.startsWith(">")) {
					String id = line.substring(1);
					ids.add(id);
					if (checkBothStrands) {
						ids.add(id+"|revcom");
					}
				}
				line = line.toUpperCase();
			}
		} while (inputStream.hasNextLine() && !line.startsWith(">"));
		String thisSeq = sb.toString();
		seqs.add(thisSeq);
		if (checkBothStrands) {
			sb = sb.reverse();
			for (int i = 0; i < sb.length(); i++) {
				switch (sb.charAt(i)) {
				case '0': sb.setCharAt(i, '3'); break;
				case '1': sb.setCharAt(i, '2'); break;
				case '2': sb.setCharAt(i, '1'); break;
				case '3': sb.setCharAt(i, '0'); break;
				}	
			}
			seqs.add(sb.toString());
		}
	}

	return new FastaSequences(ids, seqs);
}



public static MotifInstance scanMotifMax(String seq, double[][] logScoringMatrix) {
	int maxPosition = -1;
	double maxScore = -9999;
	int motifWidth = logScoringMatrix[0].length;

	if (seq.contains("A") || seq.contains("G") || seq.contains("C") || seq.contains("T")) {
		seq = SeqTools.convertToDigitRepresentation(seq);
	}

	int validWindowCounter = 0;

	for (int i = 0; i < seq.length() - motifWidth + 1; i++) {
		String kmer = seq.substring(i, i + motifWidth);
		if (kmer.contains("N")) {
			continue;
		}
		validWindowCounter++;
		char[] chars = kmer.toCharArray();
		double thisScore = 0;
		for (int j = 0; j < chars.length; j++) {
			thisScore += logScoringMatrix[chars[j]-48][j];
		} 
		if (thisScore > maxScore) {
			maxScore = thisScore;
			maxPosition = i;
		}
	}
	if ((maxPosition == -1) || (validWindowCounter < 5)) {
		return null;
	} else {
		return new MotifInstance("", "", maxPosition, maxScore); 
	}
}

public static ArrayList<MotifInstance> scanMotif(String id, String seq, double[][] logScoringMatrix, double threshold, ArrayList<String> siteset) {
	ArrayList<MotifInstance> instances = new ArrayList<MotifInstance>();
	int motifWidth = logScoringMatrix[0].length;
	if (seq.contains("A") || seq.contains("G") || seq.contains("C") || seq.contains("T")) {
		seq = convertToDigitRepresentation(seq);
	}
	for (int i = 0; i < seq.length() - motifWidth + 1; i++) {
		String kMer = seq.substring(i, i+motifWidth);
		if (kMer.contains("N")) {
			continue;
		}

		char[] chars = kMer.toCharArray();
		double thisScore = 0.0;
		for (int j = 0; j < chars.length; j++) {
			thisScore += logScoringMatrix[chars[j]-48][j];
		}
		if (thisScore > threshold) {
			instances.add(new MotifInstance(id, convertToDNARepresentation(kMer), i, thisScore));
		}
	}
	return instances;
}
public static double getMaxScoreForALogScoringMatrix(double[][] logScoringMatrix) {
	int k = logScoringMatrix[0].length;
	double maxScore = 0.0;
	for (int i = 0; i < k; i++) {
		double thisMax = -999999.0;
		for (int j = 0; j < 4; j++) {
			if (logScoringMatrix[j][i] > thisMax) {
				thisMax = logScoringMatrix[j][i];
			}
		}
		maxScore += thisMax;
	}
	return maxScore;
}

public static ArrayList<MotifInstance> findMotifInstances(Pwm theta, Pwm background, ArrayList<String> ids, ArrayList<String> seqs, ArrayList<String> siteset) {
	ArrayList<MotifInstance> results = new ArrayList<MotifInstance>();
	double[][] logScoringMatrix = getLogScoringMatrix(theta, background);
	double maxPossibleScore = getMaxScoreForALogScoringMatrix(logScoringMatrix);
	for (int i = 0; i < ids.size(); i++) {
		String id = ids.get(i);
		String seq = seqs.get(i);
		ArrayList<MotifInstance> instancesOnOneSeq = scanMotif(id, seq, logScoringMatrix, maxPossibleScore * SCORING_CUTOFF_PERCENTAGE, siteset);
		results.addAll(instancesOnOneSeq);
	}
	return results;
}

public static ArrayList<String> selectAndRemoveForTesting(ArrayList<String> seqs, double fraction) {
	ArrayList<String> sequencesSavedForTesting = new ArrayList<String>();
	//int n = (int) Math.round(seqs.size() * fraction);
	//SecureRandom r = new SecureRandom();
	for (int i = 0; i < (int)Math.round(seqs.size() * fraction); i++) {
		//String s = seqs.remove(r.nextInt(seqs.size()));
		String s = seqs.remove(i);
		sequencesSavedForTesting.add(s);
	}
	return sequencesSavedForTesting;
}

public static Set<String> getSimilarKMersWithinHammingDistanceOne(String kMer) {
	Set<String> similarKMers = new HashSet<String>();
	for (int i = 0; i < kMer.length(); i++) {
		for (char j = '0'; j<='3'; j++) {
			similarKMers.add(kMer.substring(0, i) + j + kMer.substring(i+1));
		}
	}
	return similarKMers;
}

public static Set<String> getSimilarKMersShiftedByOnePosition(String kMer) {
	Set<String> similarKMersShifted = new HashSet<String>();
	for (char j = '0'; j <= '3'; j++) {
		similarKMersShifted.add(j + kMer.substring(0, kMer.length()- 1));
		similarKMersShifted.add(kMer.substring(1) + j);
	}
	return similarKMersShifted;
}

public static Set<String> getSimilarKMersWithinHammingDistanceTwo(String kMer) {
	Set<String> similarKMers = new HashSet<String>();
	for (int i = 0 ; i < kMer.length(); i++) {
		for (char changei = '0'; changei <= '3'; changei++) {
			for (int j = i+1; j < kMer.length(); j++) {
				for (char changej = '0'; changej <= '3'; changej++) {
					similarKMers.add(kMer.substring(0,i) + changei + kMer.substring(i+1,j) + changej + kMer.substring(j+1));
				}
			}	
		}

	}
	return similarKMers;
}


public static int HammingDistance(String kmer1, String kmer2) {
	int d = 0;
	for (int i = 0; i < kmer1.length(); i++) {
		if (kmer1.charAt(i)!=kmer2.charAt(i)) {
			d++;
		}
	}
	return d;
}

public static void generateDict(int l, String s, int width, ArrayList<String> dict) {
	if (l == 0) {
		dict.add(s);
		s = s.substring(0, width - 1);
	}
	else {
		String[] alphabet = {"0", "1", "2", "3"};
		for (String c : alphabet) {
			s = s + c;
			generateDict(l - 1, s, width, dict);
			s = s.substring(0, width - l);
		}
	}
}

public static Set<String> findDistantKMerSets(Set<String> kMers, int minDistance) {
	Set<String> distantKMers = new HashSet<String>();
	ArrayList<Integer> orders = new ArrayList<Integer>();
	for (int i = 0; i < kMers.size(); i++) {
		orders.add(i);
	}
	ArrayList<String> kMersArr = new ArrayList<String>(kMers);
	Collections.shuffle(orders);
	for (int i = 0; i < kMers.size(); i++) {
		if (i % 1000 == 0) {
			System.out.println(i + "/" + kMers.size());
		}
		String kMer = kMersArr.get(orders.get(i));
		boolean satisfy = true;
		for (String anotherKMer : distantKMers) {
			if (HammingDistance(kMer, anotherKMer) < minDistance) {
				satisfy = false;
				break;
			}
		}
		if (satisfy) {
			distantKMers.add(kMer);
		}

	}
	return distantKMers;
}

public static Set<String> findDistantKMerSets(ArrayList<String> kMers, int minDistance) {
	Set<String> distantKMers = new HashSet<String>();
	ArrayList<Integer> orders = new ArrayList<Integer>();
	for (int i = 0; i < kMers.size(); i++) {
		orders.add(i);
	}

	for (int i = 0; i < kMers.size(); i++) {
		if (i % 1000 == 0) {
			System.out.println(i + "/" + kMers.size());
		}
		String kMer = kMers.get(orders.get(i));
		boolean satisfy = true;
		for (String anotherKMer : distantKMers) {
			if (HammingDistance(kMer, anotherKMer) < minDistance) {
				satisfy = false;
				break;
			}
		}
		if (satisfy) {
			distantKMers.add(kMer);
		}

	}
	return distantKMers;
}

public static void generateDistantKMerList(int width, int distance) {
	//start from 11/11/4
	System.out.println("Generating dict...");
	ArrayList<String> dict = new ArrayList<String>();
	generateDict(width, "" , width, dict);
	Set<String> temp = new HashSet<String>(); 
	for (String i : dict) {
		temp.add(i);
	}
	System.out.println(dict.size());
	System.out.println("Calculting distant kmer sets...");
	Set<String> h = findDistantKMerSets(temp, distance);
	PrintWriter p = null;
	try {
		p = new PrintWriter("./offline/"+width+"mer-max"+distance+".txt");
	} catch (Exception e) {
		e.printStackTrace();
	}
	String[] distantKMerSets = new String[h.size()];
	int counter = 0;
	for (String i : h) {
		distantKMerSets[counter] = i;
		p.println(counter + "\t" + i);
		counter++;
	}
	p.close();
	System.out.println("done.Getting matched kmer list...");
	Runtime.getRuntime().gc();
	try {
		p = new PrintWriter("./offline/"+width+"mer-max"+distance+"-closest.txt");
	} catch (Exception e) {
		e.printStackTrace();
	}
	KMerEnumerator kMerEnum = new KMerEnumerator(width);
	String i;
	while ((i = kMerEnum.getNext()) != null) {
		int shortest = 9999;

		Hashtable<Integer, Integer> distances = new Hashtable<Integer, Integer>();
		for (int j = 0; j < distantKMerSets.length; j++) {
			String kMerInList = distantKMerSets[j];
			int d = HammingDistance(i,kMerInList);
			distances.put(j, d);
			if (d < shortest) {
				shortest = d;
			}
		}

		Set<Integer> mappings = new HashSet<Integer>();
		for (int j = 0; j < distantKMerSets.length; j++) {
			if (distances.get(j) == shortest) {
				mappings.add(j);
			}
		}

		for (Integer k:mappings) {
			p.print(k + "\t");
		}
		p.print("\n");
		Runtime.getRuntime().gc();
	}
	p.close();
}

public static ArrayList<String> loadPreComputedMutuallyDistantKMerList(String offlineMutuallyDistantKMerFilename) throws Exception {
	Scanner inputStream = null;
	try {
		inputStream = new Scanner(new FileInputStream(offlineMutuallyDistantKMerFilename));
	} catch (Exception err) {
		throw err;
	}
	ArrayList<String> mutuallyDistantKMers = new ArrayList<String>();
	while (inputStream.hasNextLine()) {
		String line = inputStream.nextLine().trim();
		String[] fields = line.split("\t");
		mutuallyDistantKMers.add(Integer.parseInt(fields[0]), fields[1]);
	}
	inputStream.close();
	return mutuallyDistantKMers;
}

public static ArrayList<Integer> getDifferentPositions(String str1, String str2) {
	ArrayList<Integer> diffPos = new ArrayList<Integer>();
	char[] str1chars = str1.toCharArray();
	char[] str2chars = str2.toCharArray();
	for (int i = str1.length() - 1; i >=0; i--) {
		if (str1chars[i]!=str2chars[i]) {
			diffPos.add(i);
		} else {
			return diffPos;
		}
	}
	return diffPos;

}



public static void main(String[] args)  {
	/*Pwm theta = new Pwm("rds1.pwm");
		String dir = "./data/real/p53/coordinate_convert/";
		Pwm bg = new Pwm("background.pwm");
		double pValue = 0.0;
		try {
			ArrayList<String> posSeqs = loadSeqs(dir+"only_pet_sequences.fasta", true);
			ArrayList<String> negSeqs = loadSeqs(dir+"control_sequences.fasta", true);
			pValue = rankSumTest(posSeqs, negSeqs, theta, bg);
		} catch (Exception e) {
			e.printStackTrace();
		}*/
/*	Pwm theta = new Pwm("d:/test.pwm");
	Pwm bg = new Pwm("d:/bg.pwm");
	double[][] logScoringMatrix = getLogScoringMatrix(theta, bg);
	double maxPossibleScore = getMaxScoreForALogScoringMatrix(logScoringMatrix);
	System.out.println(maxPossibleScore * SCORING_CUTOFF_PERCENTAGE);	
	System.out.println(maxPossibleScore);
*/
	FastaSequences posInputs = null;
	FastaSequences negInputs = null;
	try {
		posInputs = loadSeqs("./data/real/p53/both-only-in-neg/sequences/forward/6kq-pbmc_pos.fasta", true);
		negInputs = loadSeqs("./data/real/p53/both-only-in-neg/sequences/forward/6kq-pbmc_neg.fasta", true);
	} catch (Exception e) {
		
	}
	ArrayList<String> posSeqs = posInputs.getSeqs();
	ArrayList<String> negSeqs = negInputs.getSeqs();
	Pwm bgPos = calculateBackground(posSeqs);
	Pwm bgNeg = calculateBackground(negSeqs);
	System.out.println("Positive\n" + bgPos);
	System.out.println("Negative\n" + bgNeg);


}


}