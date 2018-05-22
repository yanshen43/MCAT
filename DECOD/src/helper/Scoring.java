package helper;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Enumeration;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Set;

public class Scoring {
	public static double diffScoreSimple(double pEst, Pwm theta, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> bgEmitProbs) {
		//For one background
		double score = 0;
		Enumeration<String> e = kMerFreqs.keys();
		while (e.hasMoreElements()) {
			String a = e.nextElement();
			double p_a = theta.emit(a);
			score += kMerFreqs.get(a)[2] * (pEst * p_a) / (pEst * p_a + (1-pEst) * bgEmitProbs.get(a)[1]);
		}
		return score;
	}


	public static double diffScoreSimple(double[] pEsts, Pwm theta, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> bgEmitProbs) {
		//For two backgrounds
		double scores[] = {0.0, 0.0};
		Enumeration<String> e = kMerFreqs.keys();
		while (e.hasMoreElements()) {
			String a = e.nextElement();
			double p_a = theta.emit(a);
			for (int i = 0; i <= 1; i++) {
				scores[i] += kMerFreqs.get(a)[i] * (pEsts[i] * p_a) / (pEsts[i] * p_a + (1 - pEsts[i]) * bgEmitProbs.get(a)[i]);
			}
		}
		return scores[1] - scores[0];
	}

	
	public static double[][] diffScoreSimpleDerivative(double pEst, Pwm theta, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> bgEmitProbs) {
		//One background
		int k = theta.getWidth();
		double[][] derivative = new double[4][k];

		Enumeration<String> e = kMerFreqs.keys();
		while (e.hasMoreElements()) {
			String a = e.nextElement();
			double theta_a = theta.emit(a);
			double b_a = bgEmitProbs.get(a)[1];
			double baseScore = (theta_a * b_a) / Math.pow((pEst * theta_a + (1 - pEst) * b_a), 2);

			for (int j=0; j<k; j++) {
				int letter = a.charAt(j) - 48;
				derivative[letter][j] += kMerFreqs.get(a)[2]* (baseScore / theta.getElement(letter, j));
			}
		}
		double header = pEst * (1-pEst);
		for (int j=0; j<k; j++) {
			for (int i=0; i<4; i++) {
				derivative[i][j] *= header;
			}
		}
		return derivative;
	}
	
	public static double[][] diffScoreSimpleDerivative(double[] pEsts, Pwm theta, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> bgEmitProbs) {
		//Two backgrounds
		int k = theta.getWidth();
		double[][][] rawDerivatives = new double[4][k][2];
		double[][] derivative = new double[4][k];

		Enumeration<String> e = kMerFreqs.keys();
		while (e.hasMoreElements()) {
			String a = e.nextElement();
			double theta_a = theta.emit(a);
			double[] b_as = bgEmitProbs.get(a);
			double[] freqs = kMerFreqs.get(a);
			double[] baseScores = {0, 0};
			for (int i = 0; i <= 1; i++) {
				baseScores[i] = (theta_a * b_as[i]) / Math.pow((pEsts[i] * theta_a + (1 - pEsts[i]) * b_as[i]), 2);
			}
			for (int j=0; j<k; j++) {
				int letter = a.charAt(j) - 48;
				for (int i = 0; i <= 1; i++) {
					rawDerivatives[letter][j][i] += freqs[i] * (baseScores[i] / theta.getElement(letter, j)); 
				}
			}
		}
		double[] headers = {pEsts[0] * (1 - pEsts[0]), pEsts[1] * (1 - pEsts[1])};
		for (int j=0; j<k; j++) {
			for (int i=0; i<4; i++) {
				derivative[i][j] = rawDerivatives[i][j][1] * headers[1] - rawDerivatives[i][j][0] * headers[0];
			}
		}
		return derivative;
	}

	//Note that the As are related to the PWMs bg and theta - be careful that this A and the PWMs must be matched!!! 
	public static Hashtable<String, Double> calculateA(Pwm bg, Pwm theta, Hashtable<String, double[]> kMerFreqs, boolean reduceCalculationSpace, int MOTIF_WIDTH) throws InterruptedException{
		//One background
		if (Thread.currentThread().isInterrupted()) {
			throw new InterruptedException();
		}
		Hashtable<String, Double> A = new Hashtable<String, Double>();
		int k = theta.getWidth();
		if (reduceCalculationSpace || MOTIF_WIDTH > 15) {
			for (String a : kMerFreqs.keySet()) {
				double X = theta.emit(a), Y = 0, Z = 0;
				for (int l = 1; l < k; l++) {
					double part2 = theta.partEmit(1, k-l, a.substring(l));
					double part4 = theta.partEmit(k-l+1, k, a.substring(0,l));
					double part1 = bg.emit(a.substring(0,l));
					double part3 = bg.emit(a.substring(l));
					Y += part1*part2;
					Z += part3*part4;

				}
				double score = X + Y + Z;
				A.put(a, score);
			}			
		}
		else {
			KMerEnumerator kMerEnum = new KMerEnumerator(k);
			String previousKMer = null, thisKMer = null;
			double[] probs = new double[2 * k - 1];
			while ((thisKMer = kMerEnum.getNext())!=null) {
				
				if (previousKMer == null) { //initialize
					previousKMer = thisKMer;
					probs[k-1] = theta.emit(thisKMer);
					for (int l = 1; l < k; l++) {
						double part2 = theta.partEmit(1, k-l, thisKMer.substring(l));
						double part4 = theta.partEmit(k-l+1, k, thisKMer.substring(0,l));
						double part1 = bg.emit(thisKMer.substring(0,l));
						double part3 = bg.emit(thisKMer.substring(l));
						probs[k-l-1] = part1*part2;
						probs[2*k-1-l] = part3*part4;
						//System.out.println(k-l-1, " ", 2*k-1-l);
					}
				} else {
					ArrayList<Integer> differentPositions = SeqTools.getDifferentPositions(thisKMer, previousKMer);
					for (int j : differentPositions) {
						int pos=0;
						for (int l = (k-1-j); l <2*k-1-j; l++) {
							probs[l] = probs[l]/theta.getElement(previousKMer.charAt(j)-48,pos) * theta.getElement(thisKMer.charAt(j)-48,pos);
							pos++;
						}
						for (int l = 0; l < k-1-j; l++) {
							probs[l] = probs[l]/bg.getElement(previousKMer.charAt(j)-48, 0) * bg.getElement(thisKMer.charAt(j)-48, 0); 
						}
						for (int l = 2*k-1-j; l<probs.length; l++) {
							probs[l] = probs[l]/bg.getElement(previousKMer.charAt(j)-48, 0) * bg.getElement(thisKMer.charAt(j)-48, 0);
						}
					}	

				}
				
				
				previousKMer = new String(thisKMer);
				if (!kMerFreqs.containsKey(thisKMer)) {
					continue;
				} 
				double score = 0;
				for (int i = 0; i < probs.length; i++) {
					score+=probs[i];
				}
				
				A.put(thisKMer, score);

			}
		}
		return A;	
	}

	
	public static Hashtable<String, double[]> calculateA(ArrayList<Pwm> bgs, Pwm theta, Hashtable<String, double[]> kMerFreqs, boolean reduceCalculationSpace, int MOTIF_WIDTH) throws InterruptedException{
		//Two backgrounds
		if (Thread.currentThread().isInterrupted()) {
			throw new InterruptedException();
		}

		Hashtable<String, double[]> A = new Hashtable<String, double[]>();
		int k = theta.getWidth();
		if (reduceCalculationSpace || MOTIF_WIDTH > 15) {
			for (String a : kMerFreqs.keySet()) {
				double X = theta.emit(a);
				double[] Y = {0.0, 0.0} , Z = {0.0, 0.0};
				for (int l = 1; l < k; l++) {
					double part2 = theta.partEmit(1, k-l, a.substring(l));
					double part4 = theta.partEmit(k-l+1, k, a.substring(0,l));
					for (int sign = 0; sign <= 1; sign++) {
						double part1 = bgs.get(sign).emit(a.substring(0,l));
						double part3 = bgs.get(sign).emit(a.substring(l));
						Y[sign] += part1*part2;
						Z[sign] += part3*part4;
					}
				}
				double[] scores = {X + Y[0] + Z[0], X + Y[1] + Z[1]};
				A.put(a, scores);

			}
		}
		else {
			KMerEnumerator kMerEnum = new KMerEnumerator(k);
			String previousKMer = null, thisKMer = null;

			double[][] probsForEachModel = new double[2][2 * k - 1];
			while ((thisKMer = kMerEnum.getNext())!=null) {
				if (previousKMer == null) { //initialize
					previousKMer = thisKMer;
					probsForEachModel[0][k-1] = theta.emit(thisKMer);
					probsForEachModel[1][k-1] = probsForEachModel[0][k-1];

					for (int l = 1; l < k; l++) {
						double part2 = theta.partEmit(1, k-l, thisKMer.substring(l));
						double part4 = theta.partEmit(k-l+1, k, thisKMer.substring(0,l));
						for (int sign = 0; sign <= 1; sign++) {
							double part1 = bgs.get(sign).emit(thisKMer.substring(0,l));
							double part3 = bgs.get(sign).emit(thisKMer.substring(l));
							probsForEachModel[sign][k-l-1] = part1*part2;
							probsForEachModel[sign][2*k-1-l] = part3*part4;
						}
					}
				} else {
					ArrayList<Integer> differentPositions = SeqTools.getDifferentPositions(thisKMer, previousKMer);
					for (int sign = 0; sign <=1; sign++) {
						Pwm thisBg = bgs.get(sign);
						for (int j : differentPositions) {
							int pos=0;
							for (int l = (k-1-j); l <2*k-1-j; l++) {
								probsForEachModel[sign][l] = probsForEachModel[sign][l]/theta.getElement(previousKMer.charAt(j)-48,pos) * theta.getElement(thisKMer.charAt(j)-48,pos);
								pos++;
							}
							for (int l = 0; l < k-1-j; l++) {
								probsForEachModel[sign][l] = probsForEachModel[sign][l]/thisBg.getElement(previousKMer.charAt(j)-48, 0) * thisBg.getElement(thisKMer.charAt(j)-48, 0); 
							}
							for (int l = 2*k-1-j; l<probsForEachModel[sign].length; l++) {
								probsForEachModel[sign][l] = probsForEachModel[sign][l]/thisBg.getElement(previousKMer.charAt(j)-48, 0) * thisBg.getElement(thisKMer.charAt(j)-48, 0);
							}
						}
					}

				}
				previousKMer = new String(thisKMer);
				if (!kMerFreqs.containsKey(thisKMer)) {
					continue;
				} 
				
				double[] scores = new double[2];
				for (int sign = 0; sign <= 1; sign++) {
					for (int i = 0; i < probsForEachModel[0].length; i++) {
						scores[sign]+=probsForEachModel[sign][i];
					}
				}
				A.put(thisKMer, scores);
				
			}
		}
		return A;
	}

	
	public static double diffScoreFull(Pwm bg, Pwm theta, double pEst, Hashtable<String, double[]> kMerFreqs, Hashtable<String, Double> A, Hashtable<String, double[]> bgEmitProbs){
		//One background
		double score = 0;
		int k = theta.getWidth();
		Set<String> kMersInEffect = new HashSet<String>();
		kMersInEffect = kMerFreqs.keySet();
		for (String a : kMersInEffect) {
			if (kMerFreqs.get(a)[1] == kMerFreqs.get(a)[0]) {
				continue;
			}
			double thisA = A.get(a);
			score += kMerFreqs.get(a)[2] * (pEst * thisA / (pEst * thisA + (1 - (2*k-1) * pEst) * bgEmitProbs.get(a)[1]));
		}
		return score;
	} 
	


	public static double diffScoreFull(ArrayList<Pwm> bgs, Pwm theta, double[] pEsts, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> As, Hashtable<String, double[]> bgEmitProbs){
		//Two backgrounds
		int k = theta.getWidth();
		double[] scores = {0, 0};
		for (String a : kMerFreqs.keySet()) {
			double[] thisAs = As.get(a);
			double[] freqs = kMerFreqs.get(a);
			double[] thisBgs = bgEmitProbs.get(a);
			for (int i = 0; i <= 1; i++) {
				double numerator = pEsts[i] * thisAs[i];
				scores[i] += freqs[i] * (numerator / (numerator + (1 - (2 * k - 1) * pEsts[i])* thisBgs[i]));
			}
			
		}
		return scores[1] - scores[0];
	} 


	public static double[][] diffScoreFullDerivative(Pwm bg, Pwm theta, double pEst, Hashtable<String, double[]> kMerFreqs, Hashtable<String, Double> A, Hashtable<String, double[]> bgEmitProbs){
		//One background
		int k = theta.getWidth();
		double[][] derivative = new double[4][k];
		Set<String> kMersInEffect = new HashSet<String>();
		kMersInEffect = kMerFreqs.keySet();
		for (String a : kMersInEffect) {
			char[] letters = a.toCharArray();
			double b_a = bgEmitProbs.get(a)[1];
			double thisA = A.get(a);
			double baseScore = kMerFreqs.get(a)[2] * b_a / Math.pow(pEst * thisA + (1 - (2*k-1) * pEst) * b_a, 2);
			for (int m = 0; m<4; m++) { //m and n are zero-based
				for (int n = 0; n<k; n++) {
					double dX = 0, dY = 0, dZ = 0;
					if (letters[n]-48 == m) {
						dX = theta.emit(a)/theta.getElement(m, n);
					}

					//dY
					for (int l = 1; l <= k - (n+1); l++) {
						if (letters[l+n]-48 == m) {
							dY += bg.emit(a.substring(0,l)) * theta.partEmit(1, k-l, a.substring(l)) / theta.getElement(m, n);
						}
					}

					//dZ
					for (int l = k - (n+1) + 1; l < k; l++) {
						if (letters[n-k+l]-48 == m) {
							dZ += bg.emit(a.substring(l)) * theta.partEmit(k-l+1, k, a.substring(0,l)) / theta.getElement(m,n);
						}
					}

					//sum
					double dA = dX + dY + dZ; //with respect to theta_mn
					derivative[m][n] += baseScore * dA;
				}
			}
		}

		double header = pEst * (1-(2*k-1)*pEst);
  		for (int m=0; m<4; m++) {
			for (int n=0; n<k; n++) {
				derivative[m][n] *= header;
			}
		}
		return derivative;
	}

	public static double[][] diffScoreFullDerivative(ArrayList<Pwm> bgs, Pwm theta, double[] pEsts, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> As, Hashtable<String, double[]> bgEmitProbs){
		//Two backgrounds
		int k = theta.getWidth();
		double[][] derivative = new double[4][k];
		double[][][] rawDerivatives = new double[4][k][2];
		for (String a : kMerFreqs.keySet()) {
			double[] freqs = kMerFreqs.get(a);
			double[] b_a = bgEmitProbs.get(a);
			double[] thisA = As.get(a);
			char[] letters = a.toCharArray();
			double[] baseScores = {0, 0};
			for (int sign = 0; sign <=1 ; sign++) {
				baseScores[sign] = freqs[sign] * b_a[sign] / Math.pow(pEsts[sign] * thisA[sign] + (1 - (2*k-1) * pEsts[sign]) * b_a[sign], 2);
			}
			for (int m = 0; m<4; m++) { //m and n are zero-based
				for (int n = 0; n<k; n++) {
					double dX = 0;
					if (letters[n]-48 == m) {
						dX = theta.emit(a)/theta.getElement(m, n);
					}
					double[] dY = {0, 0}, dZ = {0, 0};

					//dY
					for (int l = 1; l <= k - (n+1); l++) {
						if (letters[l+n]-48 == m) {

							double tmp = theta.partEmit(1, k-l, a.substring(l)) / theta.getElement(m, n);
							for (int sign = 0; sign <= 1; sign++) {
								dY[sign] += bgs.get(sign).emit(a.substring(0,l)) * tmp;
							}
						}
					}

					//dZ
					for (int l = k - (n+1) + 1; l < k; l++) {
						if (letters[n-k+l]-48 == m) {
							double tmp = theta.partEmit(k-l+1, k, a.substring(0,l)) / theta.getElement(m,n);
							for (int sign = 0; sign <= 1; sign++) {
								dZ[sign] += bgs.get(sign).emit(a.substring(l)) * tmp;
							}
						}
					}
					//sum
					double[] dA = {dX + dY[0] + dZ[0], dX + dY[1] + dZ[1]}; //with respect to theta_mn
					rawDerivatives[m][n][0] += baseScores[0] * dA[0];
					rawDerivatives[m][n][1] += baseScores[1] * dA[1];

				}
			}
		}

		double[] headers = {pEsts[0] * (1-(2*k-1)*pEsts[0]), pEsts[1] * (1-(2*k-1)*pEsts[1])}; 
		for (int m=0; m<4; m++) {
			for (int n=0; n<k; n++) {
				derivative[m][n] = headers[1] * rawDerivatives[m][n][1] - headers[0] * rawDerivatives[m][n][0];
			}
		}
		return derivative;
	}

	public static double evaluatePValue(ArrayList<String> posSeqs, ArrayList<String> negSeqs, Result result, ArrayList<Pwm> bgs, int nIters) {
		double p = 0.0;
		double originalScore = result.getScore();
		Pwm theta = result.getPwm();
		for (int i = 0; i < nIters; i++) {
			ArrayList<String> allSeqsTest = new ArrayList<String>();
			ArrayList<String> posSeqsTest = new ArrayList<String>();
			ArrayList<String> negSeqsTest = new ArrayList<String>();
			for (String s : posSeqs) {
				allSeqsTest.add(s);
			}
			for (String s: negSeqs) {
				allSeqsTest.add(s);
			}
			Collections.shuffle(allSeqsTest);
			for (int j = 0; j < posSeqs.size(); j++) {
				posSeqsTest.add(allSeqsTest.get(j));
			}
			for (int j = 0; j < negSeqs.size(); j++) {
				negSeqsTest.add(allSeqsTest.get(j + posSeqs.size()));
			}
			
		}
		return p;
	}
}
