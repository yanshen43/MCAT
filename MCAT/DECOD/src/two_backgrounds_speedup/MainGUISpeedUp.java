package two_backgrounds_speedup;
import helper.FastaSequences;
import helper.MotifInstance;
import helper.Pwm;
import helper.Result;
import helper.Scoring;
import helper.SeqTools;

import java.awt.BorderLayout;
import java.awt.Container;
import java.awt.Dimension;
import java.awt.EventQueue;
import java.awt.Font;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.net.URLEncoder;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Enumeration;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.swing.BorderFactory;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import javax.swing.filechooser.FileNameExtensionFilter;

import org.apache.commons.math.stat.descriptive.SummaryStatistics;



/**
 * @author Shan
 *
 * Change logs:
 * 1.01: Released 2011.10.24
 *       - Fixed wrong formatting of decimal points in international regional settings
 *       - Allowed parameters to be set when calling the program in command line but asking to run in GUI mode
 */
public class MainGUISpeedUp {

	/**
	 * @param args
	 *  
	 */
	public static boolean calculateTheoreticalScoreOnly = false;
	public static boolean commandLineMode = false; 
	
	public static boolean reduceCalculationSpace = true;
	public static boolean reduceSearchSpace = true;

	public static boolean ignoreSimpleRepeats = true;
	//private static boolean limitKMerCountPerSequenceTo2 = true;
	private static boolean limitKMerCountPerSequence = true;
	private static Integer UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ = 2;
	
	private static boolean CHECK_BOTH_STRANDS = true;
	public static boolean removeReverseComplementOfAFoundMotif = true;
	public static boolean evaluateMotifSignificance = false;
	public static Integer N_PVALUE_ITERS = 100;
	
	public static boolean USE_TWO_BACKGROUND_MODELS = true;
	public static boolean USE_TWO_PESTS = true;
	public static final String VERSION_NUMBER = "V1.01 (2011.10.24)";
	
	
	//public static int TOP_KMER_DIFFERENCE_THRESHOLD = (int)Math.pow(4,6);//TODO
	//public static Set<String> topKMers;
	public static double pEstPreassigned = 1.0;
	public static final char[] alphabet = {'0', '1', '2', '3'}; //represents A, C, G, T internally
	public static final String STAMP_URL = "http://www.benoslab.pitt.edu/stamp/run.php";
	public static final String[] possibleDatabases = {"JASPAR v2010", "JASPAR v3 (Families labeled)", "TRANSFAC v11.3", "TRANSFAC (Families labeled)", "FlyReg (Bergman/Pollard)", "Fly (curated by Bergman)", "Arabidopsis (AGRIS)", "Plants (AthaMap)", "Plants (PLACE)", "Uniprobe (Bulyk)", "Yeast (Harbison et al.)", "Yeast (Maclsacc et al.)", "E.Coli (DPInteract)", "Prokaryotic (RegTransBase)", "All the above", "Selected Eukaryotic", "Predicted Human (Xie et al.)", "Predicted Fly (Tiffin)", "11 JASPAR v3 FBPs (Sandelin)", "17 JASPAR FBPs (Mahony)"};
	public static final String[] possibleDatabaseCodes = {"JASPAR",  "JASPAR_Fams", "TRANSFAC", "TRANSFAC_Fams", "FlyReg", "BergmanFly", "AGRIS", "AthaMap", "Place", "Uniprobe", "HarbisonYeast", "MacIsaacYeast", "DPInteract", "RegTransBase", "ALL", "Animal", "Xie", "Tiffin", "JASPAR_FBP", "Mahony_FBP"};
	private static Integer N_MOTIFS = 10;
	private static Integer MOTIF_WIDTH = 8;
	private static Integer MOTIF_CARDINALITY = 20;
	private static Integer ITERATION = 50;
	private static ArrayList<String> siteSet;
	//public static ArrayList<String> dict = new ArrayList<String>(); //stores all possible k-mers
	public static Set<String> dict = new HashSet<String>(); //stores all possible k-mers
	private static String posSeqFilename = "";
	private static String negSeqFilename = "";
	public static String trueThetaFile = "";
	private static boolean outputToFile = true;
	private static String OUTPUT_FILENAME;
	private static Pwm bgNeg = null;
	private static Pwm bgPos = null;
	private static String helpMessage = "Usage: java -jar DECOD-20111024.jar -nogui -pos <positive_sequence_file> -neg <negative_sequence_file> -w <width> -nmotif <#motifs to find> -niter <#iterations> -o <output_file> [-exact] [-withrepeat] [-strand forward|both] [-pmotif <#Motif per positive sequence>] [-countupperlimit <Upper limit of kmer count in a sequence>]\nSee README.txt for details\n";
	//private static boolean useComplicatedBackgroundModel = false;
	private static boolean verboseMode = false;
	
	private static boolean useSimpleDerivative = false;
	
	//private static boolean usePetersSuggestionForRemovingFirstMotif = true;
	private static int restrictToTopKMersWhenReEstimatingP = 100;

	private static JLabel statusMessageLabel = new JLabel();
	private static JTextArea outputTextArea = new JTextArea(6,42);
	private static JButton generateReportButton = new JButton("Generate report...");
	private static JButton sendToSTAMPButton = new JButton("Match with known motifs using STAMP...");
	private static StringBuilder outputStringBuilder = new StringBuilder();

	private static int currentMotifNumber = 1;
	private static ArrayList<Result> finalResults;

	public static void updateOutput(final String s) {
		outputStringBuilder.append(s+"\n");
		if (!commandLineMode) {
		EventQueue.invokeLater(new Runnable() {
			public void run() {
				outputTextArea.append(s+"\n");
				outputTextArea.setCaretPosition(outputTextArea.getText().length()-1);

			}
		});
		}
		else {
			System.out.println(s);
		}

	}

	private static void updateStatus(final String s) {
		if (!commandLineMode) {
			EventQueue.invokeLater(new Runnable() {
				public void run() {
					statusMessageLabel.setText(s);
				}
			});
		}
		else {
			System.out.println(s);
		}

	}
	private static boolean isDouble(String s) {
		try {
			Double.parseDouble(s);
			return true;
		} catch (Exception e) {
			return false;
		}
	}

	private static double[][] matrixMinus(Pwm pwmA, Pwm pwmB) {
		double[][] matA = pwmA.getMatrix();
		double[][] matB = pwmB.getMatrix(); 
		int m = matA.length;
		int n = matA[0].length;
		double[][] result = new double[m][n];
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				result[i][j] = matA[i][j] - matB[i][j];
			}
		}
		return result;
	}


	private static double matrixDotProduct(double[][] matA, double[][] matB) {
		int m = matA.length;
		int n = matA[0].length;
		double result = 0;
		for (int i=0; i<m; i++) {
			for (int j=0; j<n; j++) {
				result += matA[i][j] * matB[i][j];
			}
		}
		return result;
	}	
	

	private static ArrayList<String> initialize(Set<String> siteSetSearchSpace) {
		ArrayList<String> seeds = new ArrayList<String>();
		ArrayList<String> kMers = new ArrayList<String>(siteSetSearchSpace);
		SecureRandom r = new SecureRandom();
		for (int i = 0; i < MOTIF_CARDINALITY; i++) {
			boolean repeat = true;
			int index = 0;
			while (repeat) {
				index = r.nextInt(kMers.size());
				repeat = seeds.contains(kMers.get(index));
			} 
			seeds.add(kMers.get(index));
		}
		return seeds;
	}


	/*private static Result updateGradientDescent(Pwm bg, double pEst, Hashtable<String,Double> kMerFreqDiffs, Hashtable<String, Double> bgEmitProbs, Hashtable<String, Integer> posKMerCounts) {
		boolean exitIteration = false;
		Pwm theta = new Pwm(siteSet);
		Hashtable<String, Double> A_Theta = new Hashtable<String, Double>();
		double alpha = 1;
		double previousScore = -9999;
		Pwm previousTheta = null;
		while (!exitIteration) {
			A_Theta = calculateA(bg, theta, kMerFreqDiffs);
			double score = diffScoreFull(bg, theta, pEst, kMerFreqDiffs, A_Theta, bgEmitProbs);
			if (score - previousScore < 1e-10) {
				System.out.println(score);
				exitIteration = true;
			}
			else {
				previousScore = score;
				previousTheta = theta;
				System.out.println("Score="+score);
				double[][] deriv = diffScoreFullDerivative(bg, theta, pEst, kMerFreqDiffs, A_Theta, bgEmitProbs);
				double[][] derivTimesAlpha = matrixScalarProduct(deriv, alpha);
				double[][] updatedMatrix = matrixAdd(derivTimesAlpha, theta.getMatrix());
				double[][] rescaledMatrix = rescale(updatedMatrix);
				theta.setMatrix(rescaledMatrix);
			}
			//System.out.println("Theta=");
			//System.out.println(theta);
		}
		Result thisResult = new Result(previousScore, previousTheta);
		return thisResult;
	}*/

	private static Result update(ArrayList<Pwm> bgs, double[] pEsts, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> bgEmitProbs, Set<String> siteSetSearchSpace, boolean doExactCalculation) throws InterruptedException {
		boolean exitIteration = false;
		Pwm theta = new Pwm(siteSet);
		double diffScoreBeforeDel = 0;
		Hashtable<String, Double> A_Theta = new Hashtable<String, Double>();
		Hashtable<String, double[]> A_Thetas = new Hashtable<String, double[]>();
		if (useSimpleDerivative) {
			if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
				diffScoreBeforeDel = Scoring.diffScoreSimple(pEsts[1], theta, kMerFreqs, bgEmitProbs);
			} else {
				diffScoreBeforeDel = Scoring.diffScoreSimple(pEsts, theta, kMerFreqs, bgEmitProbs);
			}
		}
		else {
			if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
				A_Theta = Scoring.calculateA(bgs.get(1), theta, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
				diffScoreBeforeDel = Scoring.diffScoreFull(bgs.get(1), theta, pEsts[1], kMerFreqs, A_Theta, bgEmitProbs);
			} else {
				A_Thetas = Scoring.calculateA(bgs, theta, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
				diffScoreBeforeDel = Scoring.diffScoreFull(bgs, theta, pEsts, kMerFreqs, A_Thetas, bgEmitProbs);
			}
		}
		double previousImprovement = 9999;
		while (!exitIteration) {
			if (verboseMode) {
				updateOutput("score = "+diffScoreBeforeDel);
			}
			//Deletion
			int indexToBeRemoved = 0;
			double highestScore = -9999;
			Hashtable<String, Double> A_DelHighest = new Hashtable<String, Double>();
			Hashtable<String, double[]> A_DelHighests = new Hashtable<String, double[]>();

			if (Thread.currentThread().isInterrupted()) {
				throw new InterruptedException();
			}
			Pwm thetaBeforeDel = new Pwm(siteSet);
			double[][] deriv = new double[4][thetaBeforeDel.getWidth()];
			if (useSimpleDerivative) {
				if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
					deriv = Scoring.diffScoreSimpleDerivative(pEsts[1], thetaBeforeDel, kMerFreqs, bgEmitProbs);
				} else {
					deriv = Scoring.diffScoreSimpleDerivative(pEsts, thetaBeforeDel, kMerFreqs, bgEmitProbs);
				}
			} else {
				
				if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
					Hashtable<String, Double> A_BeforeDel = Scoring.calculateA(bgs.get(1), thetaBeforeDel, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
					deriv = Scoring.diffScoreFullDerivative(bgs.get(1), thetaBeforeDel, pEsts[1], kMerFreqs, A_BeforeDel, bgEmitProbs);
				} else {
					Hashtable<String, double[]> A_BeforeDels = Scoring.calculateA(bgs, thetaBeforeDel, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
					deriv = Scoring.diffScoreFullDerivative(bgs, thetaBeforeDel, pEsts, kMerFreqs, A_BeforeDels, bgEmitProbs);
					
				}

			}
			//estimate the difference brought about by deleting each k-mer.
			ArrayList<String> siteSetDel;
			for (int i = 0; i < siteSet.size(); i++) {
				siteSetDel = new ArrayList<String>(siteSet);
				siteSetDel.remove(i);
				Pwm thetaDel = new Pwm(siteSetDel);
				double[][] pwmDiff = matrixMinus(thetaDel, thetaBeforeDel);
				double scoreEst = matrixDotProduct(deriv, pwmDiff);
				if (scoreEst > highestScore) {
					highestScore = scoreEst;
					indexToBeRemoved = i;
				}
			}

			String deleted = siteSet.remove(indexToBeRemoved);

			//Add
			if (Thread.currentThread().isInterrupted()) {
				throw new InterruptedException();
			}
			Pwm thetaDel = new Pwm(siteSet);
			Map<String, Double> allPosKMerImprovements = new HashMap<String, Double>();

			if (useSimpleDerivative) {
				if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
					deriv = Scoring.diffScoreSimpleDerivative(pEsts[1], thetaDel, kMerFreqs, bgEmitProbs);
				} else {
					deriv = Scoring.diffScoreSimpleDerivative(pEsts, thetaDel, kMerFreqs, bgEmitProbs);
				}
			} else {
				if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)){
					A_DelHighest = Scoring.calculateA(bgs.get(1), new Pwm(siteSet), kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
					deriv = Scoring.diffScoreFullDerivative(bgs.get(1), thetaDel, pEsts[1], kMerFreqs, A_DelHighest, bgEmitProbs);
				} else {
					A_DelHighests = Scoring.calculateA(bgs, new Pwm(siteSet), kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
					deriv = Scoring.diffScoreFullDerivative(bgs, thetaDel, pEsts, kMerFreqs, A_DelHighests, bgEmitProbs);
					
				}
			}


			//sort the hash table by values

			//new, do not do exact calculation. only use estimations.
			if (doExactCalculation) {
				for (String a : siteSetSearchSpace) {
					ArrayList<String> siteSetAdd = new ArrayList<String>(siteSet);
					siteSetAdd.add(a);
					Pwm thetaAdd = new Pwm(siteSetAdd);
					double[][] pwmDiff = matrixMinus(thetaAdd, thetaDel);
					double diffEst = matrixDotProduct(deriv, pwmDiff);
					if (diffEst >= -1e-5) {
						allPosKMerImprovements.put(a, diffEst);
					}
				}

				List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(allPosKMerImprovements.entrySet());
				Collections.sort(list, new Comparer());
				int counter = 0;
				exitIteration = true;
				for (Map.Entry<String, Double> e: list) {
					counter++;
					String s = e.getKey();
					ArrayList<String> siteSetAdd = new ArrayList<String>(siteSet);
					siteSetAdd.add(s);
					Pwm thetaAdd = new Pwm(siteSetAdd);
					double diffExact = 0;

					if (useSimpleDerivative) {
						if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
							diffExact = Scoring.diffScoreSimple(pEsts[1], thetaAdd, kMerFreqs, bgEmitProbs);
						} else {
							diffExact = Scoring.diffScoreSimple(pEsts, thetaAdd, kMerFreqs, bgEmitProbs);
						}
					} else {
						if ((!USE_TWO_BACKGROUND_MODELS) && (!USE_TWO_PESTS)) {
							Hashtable<String, Double> A_Add = Scoring.calculateA(bgs.get(1), thetaAdd, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
							diffExact = Scoring.diffScoreFull(bgs.get(1), thetaAdd, pEsts[1], kMerFreqs, A_Add, bgEmitProbs);
						} else {
							Hashtable<String, double[]> A_Adds = Scoring.calculateA(bgs, thetaAdd, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
							diffExact = Scoring.diffScoreFull(bgs, thetaAdd, pEsts, kMerFreqs, A_Adds, bgEmitProbs);
						}
					}

					if (diffExact > diffScoreBeforeDel) {
						siteSet.add(s);
						diffScoreBeforeDel = diffExact;
						exitIteration = false;
						break;
					}
					if (counter>500) {
						break;
					}
				}
			}
			else { // do not do exact calculation
				double maxImprovement = -9999;
				String kMerWithMaxImprovement = null;
				for (String a : siteSetSearchSpace) {
					ArrayList<String> siteSetAdd = new ArrayList<String>(siteSet);
					siteSetAdd.add(a);
					Pwm thetaAdd = new Pwm(siteSetAdd);
					double[][] pwmDiff = matrixMinus(thetaAdd, thetaDel);
					double diffEst = matrixDotProduct(deriv, pwmDiff);
					if (diffEst > maxImprovement) {
						maxImprovement = diffEst;
						kMerWithMaxImprovement = a;
					}
				}

				if (maxImprovement > 0 && maxImprovement!=previousImprovement) {
					siteSet.add(kMerWithMaxImprovement);
					diffScoreBeforeDel = diffScoreBeforeDel + maxImprovement;
					exitIteration = false;
					previousImprovement = maxImprovement;
				} else {
					exitIteration = true;
				}
			}
			/*			//original, do exact calculation after the estimations
  			List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(allPosKMerImprovements.entrySet());
			Collections.sort(list, new Comparer());
			int counter = 0;
			exitIteration = true;
			for (Map.Entry<String, Double> e: list) {
				counter++;
				String s = e.getKey();
				ArrayList<String> siteSetAdd = new ArrayList<String>(siteSet);
				siteSetAdd.add(s);
				Pwm thetaAdd = new Pwm(siteSetAdd);
				double diffExact = 0;

				if (useSimpleDerivative) {
					diffExact = diffScoreSimple(pEst, thetaAdd, kMerFreqs, bgEmitProbs);
				} else {
					Hashtable<String, Double> A_Add = calculateA(bg, thetaAdd, kMerFreqs);
					diffExact = diffScoreFull(bg, thetaAdd, pEst, kMerFreqs, A_Add, bgEmitProbs);
				}

				if (diffExact > diffScoreBeforeDel) {
					siteSet.add(s);
					diffScoreBeforeDel = diffExact;
					exitIteration = false;
					break;
				}
				if (counter>500) {
					break;
				}
			}*/
			if (exitIteration) {
				siteSet.add(deleted);
			}
		}

		return new Result(diffScoreBeforeDel, siteSet, new Pwm(siteSet), pEsts);
	}

	/*	private static void outputDouble2DArray(double[][] arr) {
		for (int i=0; i<arr.length; i++) {
			for (int j=0; j<arr[0].length; j++) {
				System.out.print(arr[i][j]+" ");
			}
			System.out.println();
		}
	}*/



	public static Hashtable<String, Double> scoreAllKMers(Pwm theta, Pwm background, Hashtable<String, Double> posKMerFreqs) {
		Hashtable<String, Double> scores = new Hashtable<String, Double>();
		Enumeration<String> e = posKMerFreqs.keys();
		theta.convertToScoringMatrix(background);
		while (e.hasMoreElements()) {
			String kMer = e.nextElement();
			scores.put(kMer, theta.score(kMer));
		}
		return scores;
	}

	public static Hashtable<String, Double> scoreAllKMers(Pwm theta, Hashtable<String, Double> posKMerFreqs) {
		Hashtable<String, Double> scores = new Hashtable<String, Double>();
		Enumeration<String> e = posKMerFreqs.keys();
		while (e.hasMoreElements()) {
			String kMer = e.nextElement();
			scores.put(kMer, theta.emit(kMer));
		}
		return scores;
	}



	private static double L2Norm(Hashtable<String, Double> hash) {
		Enumeration<String> e = hash.keys();
		double norm = 0;
		while (e.hasMoreElements()) {
			String kMer = e.nextElement();
			Double value = hash.get(kMer);
			norm += Math.pow(value,2);
		}
		return Math.sqrt(norm);
	}

	/*
	 * Estimate the expected frequency of the kMer including the convolution
	 */
	public static double estimateExpectedFreq(Pwm theta, Pwm bg, String kMer) {
		double X = theta.emit(kMer), Y = 0, Z = 0;
		int k = kMer.length();
		for (int l = 1; l < k; l++) {
			double part1 = bg.emit(kMer.substring(0,l));
			double part2 = theta.partEmit(1, k-l, kMer.substring(l));
			double part3 = bg.emit(kMer.substring(l));
			double part4 = theta.partEmit(k-l+1, k, kMer.substring(0,l));
			Y += part1*part2;
			Z += part3*part4;
		}
		return X + Y + Z;
	}

	public static double reEstimateP(Hashtable<String, double[]> kMerFreqs, Pwm theta, Pwm bg, Hashtable<String, Double> expectedFreqsNormalized) {

		Hashtable<String, Double> part1 = new Hashtable<String, Double>();
		Hashtable<String, Double> part2 = new Hashtable<String, Double>();

		List<Map.Entry<String, Double>> list = new ArrayList<Map.Entry<String, Double>>(expectedFreqsNormalized.entrySet());
		Collections.sort(list, new Comparer());
		int topCounter = 0;
		//Set<String> topExpectedKMers = new HashSet<String>();
		for (Map.Entry<String, Double> entry : list) {
			topCounter++;
			if (topCounter > restrictToTopKMersWhenReEstimatingP) {
				break;
			}
			String kMer = entry.getKey();
			Double expectedNormalized = entry.getValue();


			part1.put(kMer, kMerFreqs.get(kMer)[1] - kMerFreqs.get(kMer)[0]);
			part2.put(kMer, expectedNormalized - kMerFreqs.get(kMer)[0]);
		}

		return L2Norm(part1)/L2Norm(part2);
	}

	public static Result run(ArrayList<Pwm> bgs, double[] pEsts, Hashtable<String, double[]> kMerFreqs, Hashtable<String, double[]> bgEmitProbs, ArrayList<String> posSeqs, ArrayList<String> negSeqs, Set<String> initialSiteSetSearchSpace, Hashtable<String, double[]> originalKMerFreqs, boolean showMessage) throws InterruptedException {
		//for one bg scenario
		Result finalResult = new Result(-9999, new ArrayList<String>(), new double[2]);
		long timeStart = System.currentTimeMillis();
		for (int r = 0; r < ITERATION; r++) {
			if (showMessage) {
				updateStatus("Motif " + currentMotifNumber + " / " + N_MOTIFS + ", iteration " + (r+1) + " / " + ITERATION );
			}
			if (verboseMode) {
				updateOutput("Iteration = "+r);
			}
			Result thisRun = null;
			if (!reduceSearchSpace) {
				siteSet = initialize(initialSiteSetSearchSpace);
				thisRun =  update(bgs, pEsts, kMerFreqs, bgEmitProbs, initialSiteSetSearchSpace, true);
			} else {
				//Crude search
				if (showMessage) {
					updateStatus("Motif " + currentMotifNumber + " / " + N_MOTIFS + ", iteration " + (r+1) + " / " + ITERATION + " Crude search...");
				}
				siteSet = initialize(initialSiteSetSearchSpace);
				Result thisRun1 = update(bgs, pEsts, kMerFreqs, bgEmitProbs, initialSiteSetSearchSpace, false);

				//Refined search
				if (showMessage) {
					updateStatus("Motif " + currentMotifNumber + " / " + N_MOTIFS + ", iteration " + (r+1) + " / " + ITERATION + " Refined search...");
				}
				Set<String> refinedSiteSetSearchSpace = new HashSet<String>();
				Pwm currentPwm = thisRun1.getPwm();

				for (String kMer: originalKMerFreqs.keySet()) {
					if (originalKMerFreqs.get(kMer)[1] <= originalKMerFreqs.get(kMer)[0]) {
						continue;
					}
					double score = currentPwm.emit(kMer);
					if (score >= Math.pow(0.5, Math.ceil(MOTIF_WIDTH/2)) * Math.pow(0.1, MOTIF_WIDTH - Math.ceil(MOTIF_WIDTH/2))) {
						refinedSiteSetSearchSpace.add (kMer);
					}
				}
	
				if (refinedSiteSetSearchSpace.size() <= MOTIF_CARDINALITY) {
					thisRun = thisRun1;
				} else {
					siteSet = initialize(refinedSiteSetSearchSpace);
					Result thisRun2 = update(bgs, pEsts, kMerFreqs, bgEmitProbs, refinedSiteSetSearchSpace, true);
					if (thisRun1.getScore() > thisRun2.getScore()) {
						thisRun = thisRun1;
					} else {
						thisRun = thisRun2;
					}
				}
			}
			if (verboseMode) {
				updateOutput("Recovered PWM:");
				updateOutput(thisRun.getPwm().toString());
			}

			if ((thisRun.getScore()) > finalResult.getScore()) {
				finalResult = thisRun;
				
			}
		}

		/*		if (reEstimatePEst) {
			updateOutput("Result before re-estimating P:");
			updateOutput("Score = "+finalResult.getScore());
			updateOutput("PWM =");
			updateOutput(finalResult.getPwm().toString());
			//re-estimate pEst and run one more iteration with the new pEst
			double pEstNew = reEstimateP(kMerFreqs, finalResult.getPwm(), bg);
			updateOutput("old pEst = " + pEst + ", new pEst = " + pEstNew);
			siteSet = finalResult.getSites();
			finalResult = update(bg, pEstNew, kMerFreqs, bgEmitProbs, null);
		}

		 */		
		long timeTake = System.currentTimeMillis()-timeStart;

		updateOutput(finalResult.getPwm().toString());
		updateOutput("#Score = "+finalResult.getScore());
		updateOutput("#Time = " + timeTake/1000.0+ "seconds");
		
		if (evaluateMotifSignificance) {
			/*Double pValue = Double.NaN;
			try {
				if (!USE_TWO_BACKGROUND_MODELS) {
					pValue = SeqTools.rankSumTest(posSeqsTesting, negSeqsTesting, finalResult.getPwm(), bgs.get(0));
				} else {
					pValue = SeqTools.rankSumTest(posSeqsTesting, negSeqsTesting, finalResult.getPwm(), bgs.get(1), bgs.get(0));
				}
			} catch (Exception e) {
				e.printStackTrace();
			}

			updateOutput("#p-value = " + pValue);
			finalResult.setPValue(pValue);*/
		}	
		
		return finalResult;
	}
	

	public static double getPEst(ArrayList<String> seqs) {
		int sum = 0;
		for (String s : seqs) {
			sum += s.length();
		}
		if (CHECK_BOTH_STRANDS) {
			return (double)seqs.size()/sum/2;
		} else {
			return (double)seqs.size()/sum;
		}
	}
	


	private static Hashtable<String, Double> getNormalizedExpectedFreqs(Hashtable<String, double[]> kMerFreqs, Pwm theta, Pwm bg) {
		Hashtable<String, Double> expectedFreqs = new Hashtable<String, Double>();
		double pKMerSum = 0;
		for (String kMer : kMerFreqs.keySet()) {
			if (kMerFreqs.get(kMer)[1] == 0) {
				continue;
			}
			double pKMer = estimateExpectedFreq(theta, bg, kMer);
			pKMerSum += pKMer;
			expectedFreqs.put(kMer, pKMer);
		}
		for (String kMer : kMerFreqs.keySet()) {
			if (kMerFreqs.get(kMer)[1] == 0) {
				continue;
			}
			expectedFreqs.put(kMer, expectedFreqs.get(kMer)/pKMerSum);
		}
		return expectedFreqs;
	}

	private static void removeSignalsOfAPwm(Pwm theta, Hashtable<String, double[]> kMerFreqs, Pwm bg) {
		Hashtable<String, Double> expectedFreqsNormalized = getNormalizedExpectedFreqs(kMerFreqs, theta, bg);
		double pEstNew = reEstimateP(kMerFreqs, theta, bg, expectedFreqsNormalized);
		for (String kMer : kMerFreqs.keySet()) {
			if (kMerFreqs.get(kMer)[1] == 0) {
				continue;
			}
			double[] freqs = kMerFreqs.get(kMer);
			double expectedNormalized = expectedFreqsNormalized.get(kMer);

			freqs[1] = kMerFreqs.get(kMer)[1] - pEstNew *  expectedNormalized;
			if (freqs[1] < 0) {
				freqs[1] = 0;
			}
			kMerFreqs.put(kMer, freqs);
		}
		SeqTools.normalize(kMerFreqs);

	}


	private static Hashtable<String, double[]> getOverloadedKMerFreqs(Hashtable<String, double[]> kMerFreqs) throws InterruptedException {
		if (!reduceCalculationSpace) {
			return kMerFreqs;
		}
		Hashtable<String, double[]> overloadedKMerFreqs = new Hashtable<String, double[]>();
		SummaryStatistics stats = new SummaryStatistics();
		for (String kMer: kMerFreqs.keySet()) {
			stats.addValue(kMerFreqs.get(kMer)[2]);
		}
		double mean = stats.getMean();
		double std = stats.getStandardDeviation();
		int start = 2;
		if (MOTIF_WIDTH >=10) {
			start = 3;
		}
		for (int stdFromMean = start; stdFromMean>=0; stdFromMean--) {
			double thresholdHigh = mean + stdFromMean * std;
			double thresholdLow = mean - stdFromMean * std;
			for(String kMer: kMerFreqs.keySet()) {
				if ((kMerFreqs.get(kMer)[2] > thresholdHigh) || (kMerFreqs.get(kMer)[2] < thresholdLow)) {
					overloadedKMerFreqs.put(kMer,kMerFreqs.get(kMer));
				}
			}
			if ((MOTIF_WIDTH >= 8 && overloadedKMerFreqs.size()>1000) || ((MOTIF_WIDTH <8) && (overloadedKMerFreqs.size()>100))) {
				//System.out.println("width=" + MOTIF_WIDTH);
				//System.out.println("stdFromMean=" + stdFromMean);
				break;
			}

		}
		return overloadedKMerFreqs;

	}


	private static Set<String> getInitialSearchSpace(Hashtable<String, double[]> kMerFreqs) {
		Hashtable<String, double[]> reducedKMerFreqs = new Hashtable<String, double[]>(kMerFreqs);
		
		for (String kMer: kMerFreqs.keySet()) {
			//if ((!SeqTools.isSimpleRepeat(kMer)) && (!SeqTools.isDinucleotideRepeat(kMer))) {
				double[] freqs = kMerFreqs.get(kMer);
				reducedKMerFreqs.put(kMer, freqs);
//			}
		}
		Set<String> initialSiteSetSearchSpace = new HashSet<String>();
		
		if (!reduceSearchSpace) {
			for (String kMer : reducedKMerFreqs.keySet()) {
				if (reducedKMerFreqs.get(kMer)[2]>0) {
					initialSiteSetSearchSpace.add(kMer);
				}
			}
			return initialSiteSetSearchSpace;
		}
		//updateStatus("Finding inital reduced site set search space...");
		SummaryStatistics stats = new SummaryStatistics();
		for (String kMer: reducedKMerFreqs.keySet()) {
			stats.addValue(reducedKMerFreqs.get(kMer)[1]);
		}
		double mean = stats.getMean();
		double std = stats.getStandardDeviation();
		for (double stdFromMean = 1.0 ; stdFromMean>=0; stdFromMean = stdFromMean - 0.5) {
			double threshold = mean + stdFromMean * std;
			for(String kMer: reducedKMerFreqs.keySet()) {
				if ((reducedKMerFreqs.get(kMer)[1] > threshold) && (reducedKMerFreqs.get(kMer)[2]>0)) {
					initialSiteSetSearchSpace.add(kMer);
				}
			}
			if (initialSiteSetSearchSpace.size() > MOTIF_CARDINALITY * 5) {
				break;
			}
		}
		return initialSiteSetSearchSpace;
	}
	public static void setupGlobalParameters(String[] args) {
		int counterForEssentialArguments = 0;
		int startPosition;
		if (args[0].equals("-nogui")) {
			startPosition = 1;
		} else {
			startPosition = 0;
		}
		for (int i = startPosition; i < args.length; i++) {
			if (args[i].equals("-pos"))  {
				counterForEssentialArguments++;
				String s = args[++i];
				posSeqFilename = s;
				File f = new File(posSeqFilename);
				if (!f.exists()) {
					System.out.println("Error: positive sequence file " + posSeqFilename + " does not exist!");
					System.exit(1);
				}
				continue;
			}			
			if (args[i].equals("-neg"))  {
				counterForEssentialArguments++;
				String s = args[++i];
				negSeqFilename = s;
				File f = new File(negSeqFilename);
				if (!f.exists()) {
					System.out.println("Error: negative sequence file " + negSeqFilename + " does not exist!");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-w"))  {
				String s = args[++i];
				try {
					MOTIF_WIDTH = Integer.parseInt(s);
				} catch (NumberFormatException e) {
					System.out.println("Error: Unrecognized value for the -w parameter!");
					System.out.println("An integer value specifying the width of the motif to search for must follow the -w parameter.");
					System.out.println("If this parameter is not provided, the value is set to 8 by default.");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-nmotif"))  {
				String s = args[++i];
				try {
					N_MOTIFS = Integer.parseInt(s);
				} catch (NumberFormatException e) {
					System.out.println("Error: Unrecognized value for the -nmotif parameter!");
					System.out.println("An integer value specifying the number of motifs to search for must follow the -nmotif parameter.");
					System.out.println("If this parameter is not provided, the value is set to 10 by default.");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-niter"))  {
				String s = args[++i];
				try {
					ITERATION = Integer.parseInt(s);
				} catch (NumberFormatException e) {
					System.out.println("Error: Unrecognized value for the -niter parameter!");
					System.out.println("An integer value specifying the number of iterations to run must follow the -niter parameter.");
					System.out.println("If this parameter is not provided, the value is set to 50 by default.");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-c"))  {
				String s = args[++i];
				try {
					MOTIF_CARDINALITY = Integer.parseInt(s);
				} catch (NumberFormatException e) {
					System.out.println("Error: Unrecognized value for the -c parameter!");
					System.out.println("An integer value specifying the motif cardinality must follow the -c parameter.");
					System.out.println("If this parameter is not provided, the value is set to 20 by default.");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-o"))  {
				String s = args[++i];
				outputToFile = true;
				OUTPUT_FILENAME = s;
				counterForEssentialArguments++;
				continue;
			}
			if (args[i].equals("-exact")) {
				reduceCalculationSpace = false;
				reduceSearchSpace = false;
				continue;
			}
			if (args[i].equals("-withrepeat")) {
				ignoreSimpleRepeats = false;
				limitKMerCountPerSequence = false;
				continue;
			}
			if (args[i].equals("-countupperlimit")) {
				limitKMerCountPerSequence = true;
				try {
					UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ = Integer.parseInt(args[++i]);
				}
				catch (NumberFormatException e) {
					System.out.println("Error: Unrecognized value for the -countupperlimit parameter!");
					System.out.println("An integer value specifying the upper limit of kmer count per sequence must follow the -countupperlimit parameter.");
					System.out.println("If this parameter is not provided, the value is set to 2 by default.");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-bgm")) {
				String s = args[++i];
				if (s.equals("one")) {
					USE_TWO_BACKGROUND_MODELS = false;
				} else if (s.equals("two")) {
					USE_TWO_BACKGROUND_MODELS = true;
				} else {
					System.out.println("Unrecognized value for the -bgm parameter!");
					System.out.println("Either \"one\" or \"two\" must follow -bgm, specifying whether to use one (estimated from the positive sequences only) or two background models");
					System.out.println("By default the program uses two background models.");
					System.exit(1);
				}
				continue;
			}
			if (args[i].equals("-pmotif")) {
				String s = args[++i];
				if (s.equals("one")) {
					USE_TWO_PESTS = false;
				} else if (s.equals("two")) {
					USE_TWO_PESTS = true;
				} else if (isDouble(s)){
					try {
						pEstPreassigned = Double.parseDouble(s);
						if (pEstPreassigned < 0) {
							throw new NumberFormatException();
						}
					} catch (NumberFormatException e) {
						System.out.println("Unrecognized value for the -pmotif parameter!");
						System.out.println("A positive float specifying the assumed number of times that the motif is present in each positive sequence must follow the -pmotif parameter");
						System.out.println("If this parameter is not provided, the value is set to 1 by default (assuming the motif occurs once per positive sequence). ");
						System.exit(1);
					}
					USE_TWO_PESTS = true;
				} else {
					System.out.println("Unrecognized value for the -pmotif parameter!");
					System.out.println("A float specifying the assumed number of times that the motif is present in each positive sequence must follow the -pmotif parameter");
					System.out.println("If this parameter is not provided, the value is set to 1 by default (assuming the motif occurs once per positive sequence). ");
					System.exit(1);
				}
				continue;
				
			}
					
		
			
			if (args[i].equals("-strand")) {
				String s = args[++i];
				if (s.equals("forward")) {
					CHECK_BOTH_STRANDS = false;
					removeReverseComplementOfAFoundMotif = false;
				} else if (s.equals("both")) {
					CHECK_BOTH_STRANDS = true;
					removeReverseComplementOfAFoundMotif = true;
				} else {
					System.out.println("Unrecognized parameters for -strand switch!");
					System.exit(1);
				}
				continue;  
			}
			
			if (args[i].equals("-truetheta")) { //for calculating score only
				String s = args[++i];
				trueThetaFile = s;
				calculateTheoreticalScoreOnly = true;
				continue;
			}
			if (args[i].equals("-bgneg")) {
				String s = args[++i];
				bgNeg = new Pwm(s);
				continue;
			}
			
			if (args[i].equals("-bgpos")) {
				String s = args[++i];
				bgPos = new Pwm(s);
				continue;
			}
			if (args[i].equals("-help")) {
				System.out.println(helpMessage);
				System.exit(1);
			}
			System.out.println("Unrecognized parameters!");
			System.out.println(helpMessage);
			System.exit(1);
		
		}
		if (commandLineMode && counterForEssentialArguments < 3) {
			System.out.println("Error: At least -pos <positive sequence filename>, -neg <negative sequence filename> and -o <output filename> must be supplied for the program to run in command-line mode!"); 
			System.out.println(helpMessage);
			System.exit(1);
		}

	}
	public static void reportInstances(Pwm theta, ArrayList<Pwm> bgs, ArrayList<String> posIds, ArrayList<String> posSeqs, ArrayList<String> negIds, ArrayList<String> negSeqs, ArrayList<String> siteset) {
		ArrayList<MotifInstance> posMotifInstances = SeqTools.findMotifInstances(theta, bgs.get(1), posIds, posSeqs, siteset);
		ArrayList<MotifInstance> negMotifInstances = SeqTools.findMotifInstances(theta, bgs.get(0), negIds, negSeqs, siteset);
		updateOutput("\n#Motif instances in positive sequences: ");
		
		for (MotifInstance i: posMotifInstances) {
			updateOutput(i.toString());
		}
		updateOutput("\n#Motif instances in negative sequences:");
		for (MotifInstance i: negMotifInstances) {
			updateOutput(i.toString());
		}
		updateOutput("\n");
	}

	public static void core(String[] args) throws InterruptedException {
		if (args.length != 0) {
			setupGlobalParameters(args);
		}
		System.out.println("Parameters set...");
		updateOutput("#Programming running parameters:");
		updateOutput("#Version = DECOD " + VERSION_NUMBER);
		updateOutput("#Positive sequence file = " + posSeqFilename);
		updateOutput("#Negative sequence file = " + negSeqFilename);
		updateOutput("#Motif width = " + MOTIF_WIDTH);
		updateOutput("#Number of motifs to search for = " + N_MOTIFS);
		updateOutput("#Motif cardinality = "+MOTIF_CARDINALITY);
		updateOutput("#Number of iteration = "+ ITERATION);
		updateOutput("#Use speedup calculation = " + reduceCalculationSpace);
		updateOutput("#Check both strands = " + CHECK_BOTH_STRANDS);
		//updateOutput("#Use two background models = " + USE_TWO_BACKGROUND_MODELS);
		//updateOutput("#Use two different p_motif estimates = " + USE_TWO_PESTS);
		updateOutput("#Ignore repeats = " + ignoreSimpleRepeats);
		updateOutput("#Upper limit on k-mer count per sequence = " + (limitKMerCountPerSequence ? UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ : limitKMerCountPerSequence));
		updateOutput("#Assumed motif occurrence frequencies (per sequence): " + pEstPreassigned);

		
		updateStatus("Loading sequences...");
		
		FastaSequences posInputs;
		FastaSequences negInputs;
		try {
			posInputs = SeqTools.loadSeqs(posSeqFilename, CHECK_BOTH_STRANDS);
		} catch (FileNotFoundException e) {
			updateStatus("Fatal error: file "+posSeqFilename+" was not found or could not be opened!");
			updateOutput(e.getStackTrace().toString());
			throw new InterruptedException();			
		}
		try {
			negInputs = SeqTools.loadSeqs(negSeqFilename, CHECK_BOTH_STRANDS);
		} catch (FileNotFoundException e) {
			updateStatus("Fatal error: file "+negSeqFilename+" was not found or could not be opened!");
			updateOutput(e.getStackTrace().toString());
			throw new InterruptedException();
		}
		ArrayList<String> posSeqs = posInputs.getSeqs();
		ArrayList<String> posIds = posInputs.getIds();
		ArrayList<String> negSeqs = negInputs.getSeqs();
		ArrayList<String> negIds = negInputs.getIds();
		
		updateStatus("Getting k-mer compositions...");
		Hashtable<String,double[]> kMerFreqs = new Hashtable<String, double[]>();
		SeqTools.getKMerComposition(posSeqs, MOTIF_WIDTH, limitKMerCountPerSequence, UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ, kMerFreqs, 1, ignoreSimpleRepeats);
		SeqTools.getKMerComposition(negSeqs, MOTIF_WIDTH, limitKMerCountPerSequence, UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ, kMerFreqs, 0, ignoreSimpleRepeats);
		
		updateStatus("Normalizing k-mer counts...");
		SeqTools.normalize(kMerFreqs);
		


		//updateStatus("Calculating k-mer frequency differences...");
		Hashtable<String, double[]> overloadedKMerFreqs = getOverloadedKMerFreqs(kMerFreqs);
		//SeqTools.normalize(overloadedKMerFreqs);

		//Using two pEsts, estimated separately depending on the positive and negative sequences
		//Using one pEsts, only estimated from positive
		double[] pEsts = {0.0, 0.0};
		pEsts[0] = getPEst(negSeqs) * pEstPreassigned;
		pEsts[1] = getPEst(posSeqs) * pEstPreassigned;
		
		
		//updateOutput("#Estimated p for positive = " + pEsts[1]);
		//updateOutput("#Estimated p for negative = " + pEsts[0]);
		//System.out.println("original kmer size="+kMerFreqs.size());
		Set<String> initialSiteSetSearchSpace = getInitialSearchSpace(kMerFreqs);
		
		//System.out.println("new search space size="+initialSiteSetSearchSpace.size());
		//System.out.println("new kmer size="+overloadedKMerFreqs.size());
		
		//precalculate the probability that all kmers are emitted by background model
		//Using one background model: all bgs estimated from positive sequences;
		//Using two background model, two bgs are estimated separately
		ArrayList<Pwm> bgs = new ArrayList<Pwm>();
		if ((bgNeg == null) && (bgPos == null)) {
			bgs.add(SeqTools.calculateBackground(negSeqs));
			bgs.add(SeqTools.calculateBackground(posSeqs));
		} else {
			bgs.add(bgNeg);
			bgs.add(bgPos);
		}
		
		Hashtable<String, double[]> bgEmitProbs = new Hashtable<String, double[]>();
		Enumeration<String> e = kMerFreqs.keys();
		while (e.hasMoreElements()) {
			String kMer = e.nextElement();
			
			double[] bgEmits = {bgs.get(0).emit(kMer), bgs.get(1).emit(kMer)};
			bgEmitProbs.put(kMer, bgEmits);
		}			
		


		//The following is for comparing with known motif (theta)
		if (calculateTheoreticalScoreOnly) {
			Pwm thetaTheoreticalBest = new Pwm(trueThetaFile);
			updateOutput("True theta file = " + trueThetaFile);
			updateOutput(thetaTheoreticalBest.toString());
			
			
			//A_Theta = Scoring.calculateA(bg, theta, kMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
			//diffScoreBeforeDel = Scoring.diffScoreFull(bg, theta, pEst, kMerFreqs, A_Theta, bgEmitProbs);

			double best = 0.0;
			if (!USE_TWO_BACKGROUND_MODELS) {
				Hashtable<String, Double> A = Scoring.calculateA(bgs.get(1), thetaTheoreticalBest, overloadedKMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
				best = Scoring.diffScoreFull(bgs.get(1), thetaTheoreticalBest, pEsts[1], overloadedKMerFreqs, A, bgEmitProbs);
			} else {
				Hashtable<String, double[]> As = Scoring.calculateA(bgs, thetaTheoreticalBest, overloadedKMerFreqs, reduceCalculationSpace, MOTIF_WIDTH);
				best = Scoring.diffScoreFull(bgs, thetaTheoreticalBest, pEsts, overloadedKMerFreqs, As, bgEmitProbs);
			}
			updateOutput("theoretical best score = " + best);
			System.exit(1);
		}

		updateOutput(">Motif1");
		finalResults = new ArrayList<Result>();
		currentMotifNumber = 1;

	
		Result thisResult = run(bgs, pEsts, overloadedKMerFreqs, bgEmitProbs, posSeqs, negSeqs, initialSiteSetSearchSpace, kMerFreqs, true);
		
		finalResults.add(thisResult);
		
		reportInstances(thisResult.getPwm(), bgs, posIds, posSeqs, negIds, negSeqs, thisResult.getSites());
		
		if (evaluateMotifSignificance) {
			double originalScore = thisResult.getScore();
			double pvalue = 0.0;
			for (int i = 0; i < N_PVALUE_ITERS; i++) {
				updateStatus("Evaluating motif significance iteration " + (i+1) + " / " + N_PVALUE_ITERS);
				Result shuffledResult = evaluateMotifSignificance(posSeqs, negSeqs, thisResult.getPwm());
				double newScore = shuffledResult.getScore();
				if (newScore > originalScore) {
					pvalue = pvalue + 1;
				}
			}
			pvalue = pvalue / N_PVALUE_ITERS;
			updateOutput("#pvalue = " + pvalue);
		}
		

		//Start checking the secondary motif
		for (int iter = 1; iter < N_MOTIFS; iter++) {
			updateOutput(">Motif" + (iter+1));
			//Remove signals of the first motif
			Pwm currentTheta = thisResult.getPwm();
			removeSignalsOfAPwm(currentTheta, kMerFreqs, bgs.get(1));
			if (removeReverseComplementOfAFoundMotif) {
				Pwm currentThetaRevComp = currentTheta.getReverseComplement();
				removeSignalsOfAPwm(currentThetaRevComp, kMerFreqs, bgs.get(1));
			}
			currentMotifNumber = iter+1;

			//kMerFreqDiffs = calculateKMerFreqDiffs(posKMerFreqs, negKMerFreqs);
			overloadedKMerFreqs = getOverloadedKMerFreqs(kMerFreqs);
			initialSiteSetSearchSpace = getInitialSearchSpace(kMerFreqs);

			thisResult = run(bgs, pEsts, overloadedKMerFreqs, bgEmitProbs, posSeqs, negSeqs, initialSiteSetSearchSpace, kMerFreqs, true);
			finalResults.add(thisResult);
			reportInstances(thisResult.getPwm(), bgs, posIds, posSeqs, negIds, negSeqs, thisResult.getSites());
			
			if (evaluateMotifSignificance) {
				double originalScore = thisResult.getScore();
				double pvalue = 0.0;
				for (int i = 0; i < N_PVALUE_ITERS; i++) {
					updateStatus("Evaluating motif significance iteration " + (i+1) + " / " + N_PVALUE_ITERS);
					Result shuffledResult = evaluateMotifSignificance(posSeqs, negSeqs, thisResult.getPwm());
					double newScore = shuffledResult.getScore();
					if (newScore > originalScore) {
						pvalue = pvalue + 1;
					}
				}
				pvalue = pvalue / N_PVALUE_ITERS;
				updateOutput("#pvalue = " + pvalue);
			}
		}

		if (outputToFile) {
			try {
				PrintWriter outputFile = new PrintWriter(OUTPUT_FILENAME);
				outputFile.println(outputStringBuilder.toString());
				outputFile.close();
			} catch (Exception error) {
				System.out.println("Error: output to file error");
				error.printStackTrace();
			}
		}

		EventQueue.invokeLater(new Runnable() {
			public void run() {
				statusMessageLabel.setText("Finished!");
				generateReportButton.setEnabled(true);
				sendToSTAMPButton.setEnabled(true);

			}
		});		
	}


	public static Result evaluateMotifSignificance(ArrayList<String> posSeqs, ArrayList<String> negSeqs, Pwm thetaOriginal) throws InterruptedException {
		
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
		
		Hashtable<String,double[]> kMerFreqsTest = new Hashtable<String, double[]>();
		SeqTools.getKMerComposition(posSeqsTest, MOTIF_WIDTH, limitKMerCountPerSequence, UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ, kMerFreqsTest, 1, ignoreSimpleRepeats);
		SeqTools.getKMerComposition(negSeqsTest, MOTIF_WIDTH, limitKMerCountPerSequence, UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ, kMerFreqsTest, 0, ignoreSimpleRepeats);
		SeqTools.normalize(kMerFreqsTest);
		


		//updateStatus("Calculating k-mer frequency differences...");
		Hashtable<String, double[]> overloadedKMerFreqsTest = getOverloadedKMerFreqs(kMerFreqsTest);
		
		double[] pEstsTest = {0.0, 0.0};
		pEstsTest[0] = getPEst(negSeqsTest) * pEstPreassigned;
		pEstsTest[1] = getPEst(posSeqsTest) * pEstPreassigned;
		
		
		
		ArrayList<Pwm> bgsTest = new ArrayList<Pwm>();
		bgsTest.add(SeqTools.calculateBackground(negSeqsTest));
		bgsTest.add(SeqTools.calculateBackground(posSeqsTest));
		
		Hashtable<String, double[]> bgEmitProbsTest = new Hashtable<String, double[]>();
		Enumeration<String> e = kMerFreqsTest.keys();
		while (e.hasMoreElements()) {
			String kMer = e.nextElement();
			double[] bgEmits = {bgsTest.get(0).emit(kMer), bgsTest.get(1).emit(kMer)};
			bgEmitProbsTest.put(kMer, bgEmits);
		}
		Set<String> initialSiteSetSearchSpaceTest = getInitialSearchSpace(kMerFreqsTest);
		return run(bgsTest, pEstsTest, overloadedKMerFreqsTest, bgEmitProbsTest, posSeqsTest, negSeqsTest, initialSiteSetSearchSpaceTest, kMerFreqsTest, false);

	}
	
	public static void writeHtmlReport(String dir) throws Exception {
		String plainTextFilename = dir + File.separatorChar + "decod_result.txt";
		PrintWriter plainText = new PrintWriter(plainTextFilename);
		plainText.println(outputStringBuilder.toString());
		plainText.close();

		String htmlFilename = dir + File.separatorChar + "decod_report.html";
		PrintWriter out = new PrintWriter(htmlFilename);
		StringBuilder sb = new StringBuilder();
		sb.append("<html>\n");
		sb.append("<head><title>DECOD: Deconvolved Discriminative Motif Finding Report</title></head>\n");
		sb.append("<body>\n");
		sb.append("<h1>Discriminative motif search report</h1><p>\n");
		sb.append("<h3>Parameters and settings</h3><p>");
		sb.append("<b>Program Version: </b>" + VERSION_NUMBER + "<br>\b");
		sb.append("<b>Positive sequence file: </b>" + posSeqFilename + "<br>\n");
		sb.append("<b>Negative sequence file: </b>" + negSeqFilename + "<br>\n");
		sb.append("<b>Motif width: </b>" + MOTIF_WIDTH + "<br>\n");
		sb.append("<b>Number of motifs: </b>" + N_MOTIFS + "<br>\n");
		sb.append("<b>Motif cardinality: </b>" + MOTIF_CARDINALITY + "<br>\n");
		sb.append("<b>Number of iterations run: </b>" + ITERATION + "<br>\n");
		sb.append("<b>Using speedup calculation: </b>");
		if (reduceCalculationSpace) {
			sb.append("Yes<br>\n");
		} else {
			sb.append("No<br>\n");
		}
		sb.append("<b>Check both strands of input sequences: </b>");
		if (CHECK_BOTH_STRANDS) {
			sb.append("Yes<br>\n");
		} else {
			sb.append("No<br>\n");
		}
		sb.append("<b>Ignore simple repeats: </b>");
		if (ignoreSimpleRepeats) {
			sb.append("Yes<br>\n");
		} else {
			sb.append("No<br>\n");
		}

		sb.append("<b>Upper limit on k-mer counts per sequence: </b>");
		if (limitKMerCountPerSequence) {
			sb.append(UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ  + "<br>\n");
		} else {
			sb.append("No<br>\n");
		}
		sb.append("<b>Assumed motif occurrence frequencies (per sequence): </b>" + pEstPreassigned + "<br>\n");

		sb.append("<p>\n");
		sb.append("<h3>Output files</h3><p>\n");
		sb.append("<b>Html report file (this file):</b> " + htmlFilename + "<br>\n");
		sb.append("<b>Plain text file: </b>" + plainTextFilename + "<br>\n");
		sb.append("You may query <a href=\"http://www.benoslab.pitt.edu/stamp/\">STAMP</a> to find matches to known motifs in the program interface directly, or you can submit the plain text result file to <a href=\"http://www.benoslab.pitt.edu/stamp/\">STAMP</a> manually to allow for more customizations.<p>\n");
		sb.append("<h3>Discriminative motifs found by the program</h3><p>\n");
		sb.append("<table border=\"1\">\n");
		if (!evaluateMotifSignificance) {
			sb.append("<tr><td align=\"center\">Number</td><td align=\"center\">Sequence log representation</td><td align=\"center\">Matrix representation</td><td align=\"center\">Score</td></tr>\n");
		} else {
			sb.append("<tr><td align=\"center\">Number</td><td align=\"center\">Sequence log representation</td><td align=\"center\">Matrix representation</td><td align=\"center\">Score</td><td align=\"center\">P-value</td></tr>\n");
		}
		for (int i = 0; i < finalResults.size(); i++) {
			Result r = finalResults.get(i);
			r.getPwm().convertToLogo(dir + File.separatorChar + "images" + File.separatorChar + "figmotif" + (i+1) + ".png");
			sb.append("<tr>\n");
			sb.append("\t<td> Motif " + (i+1) + "</td>\n");
			sb.append("\t<td><img src=\"./images/figmotif" + (i+1) + ".png\"</td>\n");
			sb.append("\t<td>");
			sb.append(r.getPwm().toHtml());
			sb.append("</td><td>");
			sb.append(String.format("%.3e", r.getScore()));
			sb.append("</td>");
			if (evaluateMotifSignificance) {
				sb.append("<td>");
				if (r.getPValue()==null) {
					sb.append("N/A");
				} else {
					sb.append(String.format("%.4e", r.getPValue()));
				}
				sb.append("</td>");
			}
			sb.append("</tr>\n");


		}
		sb.append("</table>");
		sb.append("</body>\n</html>");
		out.println(sb.toString());	
		out.close();
		if( !java.awt.Desktop.isDesktopSupported() ) {
			throw new Exception("Cannot open the generated html file in your default browser. ");
		}
		java.awt.Desktop desktop = java.awt.Desktop.getDesktop();
		try {
			java.net.URI uri = (new File(dir+File.separatorChar+"decod_report.html")).toURI(); 
			desktop.browse(uri);
		} catch (Exception e) {
			throw e;
		}


	}

	@SuppressWarnings("serial")
	public static void main(String[] args) throws InterruptedException {
		Locale.setDefault(Locale.US);
		class DiscMotifPanel extends JPanel {
			public DiscMotifPanel(String title) {
				GridBagLayout layout = new GridBagLayout();
				setLayout(layout);
				setBorder(BorderFactory.createTitledBorder(BorderFactory.createEtchedBorder(), title));
			}

		}
		class SequencePanel extends DiscMotifPanel {
			public SequencePanel(JFileChooser chooser) {
				super("Input sequences");
				JLabel posSeqFileLabel = new JLabel("Positive sequences (.fasta format):");
				JLabel negSeqFileLabel = new JLabel("Negative sequences (.fasta format):");
				posSeqTextField = new JTextField(18);
				posSeqTextField.setText(posSeqFilename);
				posSeqTextField.setToolTipText("Provide the input file in FASTA format containing positive sequences here.");
				negSeqTextField = new JTextField(18);
				negSeqTextField.setText(negSeqFilename);
				negSeqTextField.setToolTipText("Provide the input file in FASTA format containing negative sequences here.");
				JButton choosePosFileButton = new JButton("Choose...");
				JButton chooseNegFileButton = new JButton("Choose...");
				choosePosFileButton.setPreferredSize(new Dimension(90,20));
				chooseNegFileButton.setPreferredSize(new Dimension(90,20));
				choosePosFileButton.addActionListener(new ChoosePosFileAction(this, chooser));
				chooseNegFileButton.addActionListener(new ChooseNegFileAction(this, chooser));
				add(posSeqFileLabel,new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,2,5));
				add(negSeqFileLabel,new GBCHelper(0,1,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,2,5));
				add(posSeqTextField, new GBCHelper(1,0,1,1).setWeight(100,0).setInsets(0,0,2,0));
				add(negSeqTextField, new GBCHelper(1,1,1,1).setWeight(100,0).setInsets(0,0,2,0));
				add(choosePosFileButton, new GBCHelper(2,0,1,1).setInsets(2,5,2,5));
				add(chooseNegFileButton, new GBCHelper(2,1,1,1).setInsets(2,5,2,5));
			}


			class ChoosePosFileAction implements ActionListener {
				public ChoosePosFileAction(JPanel panel, JFileChooser chooser) {
					this.callingPanel = panel;
					this.chooser = chooser;
				}
				public void actionPerformed(ActionEvent event) {
					chooser.resetChoosableFileFilters();
					chooser.setFileFilter(new FileNameExtensionFilter("FASTA Files", "fasta"));
					int result = chooser.showDialog(callingPanel, "Select");
					if (result == JFileChooser.APPROVE_OPTION) {
						posSeqFilename = chooser.getSelectedFile().getPath();
						posSeqTextField.setText(posSeqFilename);	
					}

				}
				private JPanel callingPanel;
				private JFileChooser chooser;
			}

			class ChooseNegFileAction implements ActionListener{

				public ChooseNegFileAction(JPanel panel, JFileChooser chooser) {
					this.callingPanel = panel;
					this.chooser = chooser;
				}
				public void actionPerformed(ActionEvent event) {
					chooser.resetChoosableFileFilters();
					chooser.setFileFilter(new FileNameExtensionFilter("FASTA Files", "fasta"));
					int result = chooser.showDialog(callingPanel, "Select");
					if (result == JFileChooser.APPROVE_OPTION) {
						negSeqFilename = chooser.getSelectedFile().getPath();
						negSeqTextField.setText(negSeqFilename);
					}
				}
				private JPanel callingPanel;
				private JFileChooser chooser;
			}
			private JTextField posSeqTextField;
			private JTextField negSeqTextField;

		}
		class ParameterPanel extends DiscMotifPanel {
			public ParameterPanel(JFileChooser chooser) {
				super("Parameters");
				motifWidthLabel = new JLabel(MOTIF_WIDTH.toString());
				numberOfMotifsLabel = new JLabel(N_MOTIFS.toString());
				numberOfIterationsLabel = new JLabel(ITERATION.toString());
				motifCardinalityLabel = new JLabel(MOTIF_CARDINALITY.toString());

				checkBothStrandsCheckBox = new JCheckBox("Consider both strands of the input sequences");
				checkBothStrandsCheckBox.setSelected(CHECK_BOTH_STRANDS);
				checkBothStrandsCheckBox.addActionListener(new ActionListener() 
				{
					public void actionPerformed(ActionEvent arg0) {
						CHECK_BOTH_STRANDS = checkBothStrandsCheckBox.isSelected();
					}});
				checkBothStrandsCheckBox.setToolTipText("<html>Check this option to search for discrimnative motifs in both strands of<br> the positive sequences versus both strands of the negative sequences</html>");

				limitKMerCountsPerSequenceCheckBox = new JCheckBox("Limit k-mer counts per sequence to a maximum of 2");
				limitKMerCountsPerSequenceCheckBox.setSelected(limitKMerCountPerSequence);
				limitKMerCountsPerSequenceCheckBox.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						limitKMerCountPerSequence = limitKMerCountsPerSequenceCheckBox.isSelected();
						if (limitKMerCountPerSequence) {
							UPPER_LIMIT_OF_KMER_COUNT_PER_SEQ = 2;
						}
					}
				});
				limitKMerCountsPerSequenceCheckBox.setToolTipText("<html>Check this option to limit the k-mer counts per sequence to a maximum 2. <br>This may be useful for the algorithm to avoid being trapped by repeat elements.</html>");

				
				ignoreSimpleRepeatsCheckBox = new JCheckBox("Ignore simple one- or di-nucleotide repeats in the input sequences");
				ignoreSimpleRepeatsCheckBox.setSelected(ignoreSimpleRepeats);
				ignoreSimpleRepeatsCheckBox.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						ignoreSimpleRepeats = ignoreSimpleRepeatsCheckBox.isSelected();
					}
				});
				ignoreSimpleRepeatsCheckBox.setToolTipText("<html>Check this option to ignore simple repeats in the input sequences. <br> This may be useful if the input sequences still contain large simple <br>repeat regions after repeat-masked.</html>");
				
				speedupCheckBox = new JCheckBox("Speed up the calculations (recommended for k>=8)");
				speedupCheckBox.setSelected(reduceCalculationSpace);
				speedupCheckBox.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						reduceCalculationSpace = speedupCheckBox.isSelected();
						reduceSearchSpace = speedupCheckBox.isSelected();
					}
				});
				speedupCheckBox.setToolTipText("<html>Check this option to speed up the calculation by using only the most<br> informative k-mers(see paper for details). This is strongly <br>recommended when searching for long motifs longer than 8</html>");
				
				removeReverseComplementOfAFoundMotifCheckBox = new JCheckBox("Remove the reverse complements of found motifs");
				removeReverseComplementOfAFoundMotifCheckBox.setSelected(removeReverseComplementOfAFoundMotif);
				removeReverseComplementOfAFoundMotifCheckBox.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						removeReverseComplementOfAFoundMotif = removeReverseComplementOfAFoundMotifCheckBox.isSelected();
					}
				});
				removeReverseComplementOfAFoundMotifCheckBox.setToolTipText("<html>Check this option to let the program remove the reverse complement<br> of a found motif as well before looking for the next motif.<br>This can be helpful if the program is set to find motifs on both strands.</html> ");
				
				
				saveOutputToFileCheckBox = new JCheckBox("Save output to file:");
				saveOutputToFileCheckBox.setSelected(outputToFile);
				saveOutputToFileCheckBox.setToolTipText("<html>Check this option and select an output filename to save the outputs generated <br>by the program (as will appear in the Result section below) to a text file.</html>");
				saveOutputToFileCheckBox.addActionListener(new ActionListener() {
					public void actionPerformed(ActionEvent e) {
						outputToFile = saveOutputToFileCheckBox.isSelected();
						if (!outputToFile) {
							chooseOutputFilenameButton.setEnabled(false);
							outputFilenameTextField.setEnabled(false);
						} else {
							chooseOutputFilenameButton.setEnabled(true);
							outputFilenameTextField.setEnabled(true);
						}
					}
				});

				outputFilenameTextField = new JTextField(20);
				outputFilenameTextField.setEnabled(outputToFile);
				if (OUTPUT_FILENAME != null) {
					outputFilenameTextField.setText(OUTPUT_FILENAME);
				}
				chooseOutputFilenameButton = new JButton("Choose...");
				chooseOutputFilenameButton.addActionListener(new saveFileChooser(this, chooser));
				chooseOutputFilenameButton.setEnabled(outputToFile);
				chooseOutputFilenameButton.setPreferredSize(new Dimension(90,20));

				JSlider motifWidthSlider = addSlider(6,15,MOTIF_WIDTH,1,1);
				motifWidthSlider.addChangeListener(new motifWidthSliderStateChange());
				motifWidthSlider.setToolTipText("<html>Set the desired width of the motif to search for here. <br>Longer motifs will take longer running time.</html>");

				JSlider motifCardinalitySlider= addSlider(4,40,MOTIF_CARDINALITY,4,2);
				motifCardinalitySlider.addChangeListener(new motifCardinalitySliderStateChange());
				motifCardinalitySlider.setToolTipText("<html>Set the desired motif cardinality here. This controls the number of k-mers <br>from which the PWM of the motif is constructed. A higher cardinality results <br>in a motif with higher resolution, but also require longer running time.</html>");

				JSlider numberOfMotifSlider = addSlider(1,20,N_MOTIFS,4,1);
				numberOfMotifSlider.addChangeListener(new numberOfMotifSliderStateChange());
				numberOfMotifSlider.setToolTipText("<html>Set the desired number of discriminative motifs to search for here. <br>After each motif is found, it will be removed probabilistically from<br> and the program will re-run to find the next motif.</html>");

				JSlider numberOfIterationSlider = addSlider(5,50,ITERATION,5,5);
				numberOfIterationSlider.addChangeListener(new numberOfIterationsSliderStateChange());
				numberOfIterationSlider.setToolTipText("<html>Set the desired number of iterations to run here.<br> More iterations make the program less likely <br>to get stuck in a local optima, but requires longer<br> running time.</html>");

				JPanel motifWidthPanel = new JPanel();
				motifWidthPanel.setLayout(new GridBagLayout());
				motifWidthPanel.add(new JLabel("Motif width (default 8) : "), new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST));
				motifWidthPanel.add(motifWidthLabel, new GBCHelper(1,0,1,1).setAnchor(GBCHelper.WEST));
				motifWidthPanel.add(motifWidthSlider, new GBCHelper(0,1,2,1).setAnchor(GBCHelper.WEST));
				add(motifWidthPanel, new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,0,17));

				JPanel motifCardinalityPanel = new JPanel();
				motifCardinalityPanel.setLayout(new GridBagLayout());
				motifCardinalityPanel.add(new JLabel("Motif cardinality (default 20) : "), new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,0,0));
				motifCardinalityPanel.add(motifCardinalityLabel, new GBCHelper(1,0,1,1).setAnchor(GBCHelper.WEST));
				motifCardinalityPanel.add(motifCardinalitySlider, new GBCHelper(0,1,2,1).setAnchor(GBCHelper.WEST));
				add(motifCardinalityPanel, new GBCHelper(1,0,1,1).setInsets(0,5,0,10).setAnchor(GBCHelper.WEST));

				JPanel numberOfMotifsPanel = new JPanel();
				numberOfMotifsPanel.setLayout(new GridBagLayout());
				numberOfMotifsPanel.add(new JLabel("Number of motifs to find (default 10) : "), new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST));
				numberOfMotifsPanel.add(numberOfMotifsLabel, new GBCHelper(1,0,1,1).setAnchor(GBCHelper.WEST));
				numberOfMotifsPanel.add(numberOfMotifSlider, new GBCHelper(0,1,2,1).setAnchor(GBCHelper.WEST));
				add(numberOfMotifsPanel, new GBCHelper(0,1,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,0,17));

				JPanel numberOfIterationsPanel = new JPanel();
				numberOfIterationsPanel.setLayout(new GridBagLayout());
				numberOfIterationsPanel.add(new JLabel("Number of iterations (default 50) : "), new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,0,0));
				numberOfIterationsPanel.add(numberOfIterationsLabel, new GBCHelper(1,0,1,1).setAnchor(GBCHelper.WEST));
				numberOfIterationsPanel.add(numberOfIterationSlider, new GBCHelper(0,1,2,1).setAnchor(GBCHelper.WEST));
				add(numberOfIterationsPanel, new GBCHelper(1,1,1,1).setAnchor(GBCHelper.WEST).setInsets(0,5,0,10));

				JPanel otherParametersPanel = new JPanel();
				otherParametersPanel.setLayout(new GridBagLayout());
				otherParametersPanel.add(speedupCheckBox, new GBCHelper(0,0,3,1).setAnchor(GBCHelper.WEST));
				otherParametersPanel.add(limitKMerCountsPerSequenceCheckBox, new GBCHelper(0,1,3,1).setAnchor(GBCHelper.WEST));
				otherParametersPanel.add(ignoreSimpleRepeatsCheckBox, new GBCHelper(0,2,3,1).setAnchor(GBCHelper.WEST));
				otherParametersPanel.add(checkBothStrandsCheckBox, new GBCHelper(0,3,3,1).setAnchor(GBCHelper.WEST));
				otherParametersPanel.add(removeReverseComplementOfAFoundMotifCheckBox, new GBCHelper(0,4,3,1).setAnchor(GBCHelper.WEST));
				
				otherParametersPanel.add(saveOutputToFileCheckBox, new GBCHelper(0,5,1,1).setAnchor(GBCHelper.WEST));
				otherParametersPanel.add(outputFilenameTextField, new GBCHelper(1,5,1,1).setAnchor(GBCHelper.WEST).setInsets(1,2,0,5));
				otherParametersPanel.add(chooseOutputFilenameButton, new GBCHelper(2,5,1,1).setAnchor(GBCHelper.WEST));
				add(otherParametersPanel, new GBCHelper(0,4,4,1).setAnchor(GBCHelper.WEST).setInsets(0,0,5,0));

			}

			private JSlider addSlider(int min, int max, int init, int major, int minor) {
				JSlider s = new JSlider(min, max, init);
				s.setPaintTicks(true);
				s.setSnapToTicks(true);
				s.setPaintLabels(true);
				s.setMajorTickSpacing(major);
				s.setMinorTickSpacing(minor);
				s.setPreferredSize(new Dimension(250,50));
				return s;
			}

			public JTextField getOutputFilenameTextField() {
				return outputFilenameTextField;
			}
			private JLabel motifWidthLabel, numberOfMotifsLabel, numberOfIterationsLabel, motifCardinalityLabel;
			private JCheckBox checkBothStrandsCheckBox, limitKMerCountsPerSequenceCheckBox, saveOutputToFileCheckBox, speedupCheckBox, removeReverseComplementOfAFoundMotifCheckBox, ignoreSimpleRepeatsCheckBox;
			private JButton chooseOutputFilenameButton;



			private JTextField outputFilenameTextField = new JTextField(20);
			class motifWidthSliderStateChange implements ChangeListener {
				public void stateChanged(ChangeEvent event) {
					JSlider slider = (JSlider)event.getSource();
					MOTIF_WIDTH = slider.getValue();
					motifWidthLabel.setText(MOTIF_WIDTH.toString());
				}	
			}
			class motifCardinalitySliderStateChange implements ChangeListener {
				public void stateChanged(ChangeEvent event) {
					JSlider slider = (JSlider)event.getSource();
					MOTIF_CARDINALITY = slider.getValue();
					motifCardinalityLabel.setText(MOTIF_CARDINALITY.toString());
				}	
			}
			class numberOfMotifSliderStateChange implements ChangeListener {
				public void stateChanged(ChangeEvent event) {
					JSlider slider = (JSlider)event.getSource();
					N_MOTIFS = slider.getValue();
					numberOfMotifsLabel.setText(N_MOTIFS.toString());
				}	
			}
			class numberOfIterationsSliderStateChange implements ChangeListener {
				public void stateChanged(ChangeEvent event) {
					JSlider slider = (JSlider)event.getSource();
					ITERATION = slider.getValue();
					numberOfIterationsLabel.setText(ITERATION.toString());
				}	
			}
			class saveFileChooser implements ActionListener {
				public saveFileChooser(JPanel p, JFileChooser chooser) {
					this.parent = p;
					this.chooser = chooser;
				}
				public void actionPerformed(ActionEvent e) {
					//chooser = new JFileChooser();
					chooser.resetChoosableFileFilters();
					chooser.setFileFilter(new FileNameExtensionFilter("Text Files", "txt"));
					int result = chooser.showDialog(parent , "Save");
					if (result == JFileChooser.APPROVE_OPTION) {
						OUTPUT_FILENAME = chooser.getSelectedFile().getPath();
						outputFilenameTextField.setText(OUTPUT_FILENAME);	
					}

				} 
				private JPanel parent;
				private JFileChooser chooser;
			}

		}

		class ActionPanel extends JPanel {
			public ActionPanel(JTextField outputFilenameTextField, JFrame parent) {
				this.outputFilenameTextField = outputFilenameTextField;
				setLayout(new GridBagLayout());
				runButton = new JButton("Run");
				runButton.addActionListener(new RunAction(this));
				helpButton = new JButton("Help");
				t = new Thread(new CoreThread());
				helpButton.addActionListener(new HelpAction());
				aboutButton = new JButton("About");
				aboutButton.addActionListener(new AboutAction(parent));
				exitButton = new JButton("Exit");
				exitButton.addActionListener(new ExitAction(this));
				add(runButton, new GBCHelper(0,0,1,1).setInsets(10,10,10,10));
				//add(helpButton, new GBCHelper(1,0,1,1).setInsets(10,10,10,10));
				add(aboutButton, new GBCHelper(2,0,1,1).setInsets(10,10,10,10));
				add(exitButton, new GBCHelper(3,0,1,1).setInsets(10,10,10,10));

			}
			private JButton runButton, helpButton, aboutButton, exitButton;
			private JDialog dialog;
			//private boolean isRunning = false;
			private Thread t;
			private boolean taskFinished = false;
			private JTextField outputFilenameTextField;

			class RunAction implements ActionListener {
				public RunAction(JPanel p) {
					parent = p;
				}
				public void actionPerformed(ActionEvent event) {
					if (taskFinished) {
						t = new Thread(new CoreThread());
						taskFinished = false;
						outputStringBuilder = new StringBuilder();
						updateOutput(outputStringBuilder.toString());
					}
					if (!t.isAlive()) {
						File f = new File(posSeqFilename);
						if (!f.exists()) {
							JOptionPane.showMessageDialog(parent, "Error: the input positive sequence file you specified does not exist!", "Error", JOptionPane.ERROR_MESSAGE);
						} else {
							f = new File(negSeqFilename);
							if (!f.exists()) {
								JOptionPane.showMessageDialog(parent, "Error: the input negative sequence file you specified does not exist!", "Error", JOptionPane.ERROR_MESSAGE);
							} else {
								OUTPUT_FILENAME = outputFilenameTextField.getText();
								if (OUTPUT_FILENAME.equals("") && outputToFile) {
									JOptionPane.showMessageDialog(parent, "Error: you asked the program to save its output to a file but the output filename is blank!", "Error", JOptionPane.ERROR_MESSAGE);
								} else {
									outputTextArea.setText("");
									runButton.setText("Stop");
									sendToSTAMPButton.setEnabled(false);
									generateReportButton.setEnabled(false);
									t.start();		
								}

							}
						}
					}  else {
						t.interrupt();
					}


				}
				private JPanel parent;


			}

			class HelpAction implements ActionListener {

				@Override
				public void actionPerformed(ActionEvent e) {

				}

			}
			class AboutDialog extends JDialog {
				public AboutDialog(JFrame owner) {
					super(owner, "About DECOD", true);
					add(new JLabel("<html><h3><i>DECOD: Deconvolved discriminative motif finder</i></h3><hr><b>Please cite:</b> <p>Peter Huggins, Shan Zhong, Idit Shiff, et al. <i>DECOD: Fast and Accurate Discriminaite Motif Finding</i>, Bioinformatics, 2011, 27(17): 2361-7. <p><br><b>Bug report: </b><br>Please contact Shan Zhong at szhong@andrew.cmu.edu<br><br>See README.txt for more information.</html>"), BorderLayout.CENTER);
					JPanel panel = new JPanel();
					JButton ok = new JButton("OK");
					ok.addActionListener(new ActionListener(){
						public void actionPerformed(ActionEvent event) {
							setVisible(false);
						}
					});
					panel.add(ok);
					add(panel, BorderLayout.SOUTH);
					setSize(500,250);
				}
			}
			class AboutAction implements ActionListener {
				public AboutAction(JFrame p) {
					parent = p;
				}

				@Override
				public void actionPerformed(ActionEvent e) {
					if (dialog == null) {
						dialog = new AboutDialog(parent);
					}
					dialog.setLocationByPlatform(true);
					dialog.setVisible(true);
				}
				private JFrame parent;
			}


			class ExitAction implements ActionListener {
				public ExitAction(JPanel p) {
					parent = p;
				}
				public void actionPerformed(ActionEvent e) {
					int selection = JOptionPane.showConfirmDialog(parent, "Exit the program?", "Exit?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
					if (selection == JOptionPane.YES_OPTION) {
						System.exit(1);
					}
				}
				private JPanel parent;
			}

			class CoreThread implements Runnable {
				private volatile boolean threadDone = false;
				public void terminate() {
					threadDone = true;
				}

				public void run() {
					try {
						if (!Thread.currentThread().isInterrupted() && !threadDone) {
							core(new String[0]);
						}
					} catch (InterruptedException e) {
						statusMessageLabel.setText("Interrupted!");
						updateOutput("Interrupted!");
						Thread.currentThread().interrupt();
					} catch (Exception e) {
						statusMessageLabel.setText("Program exited with error.");
						updateOutput("Program exited with error!");
						updateOutput(e.toString());
						e.printStackTrace();
						Thread.currentThread().interrupt();
					} finally {
						this.terminate();
						runButton.setText("Run");
						taskFinished = true;
					}


				}
			}

		}

		class StatusPanel extends JPanel {
			public StatusPanel() {
				add(statusMessageLabel);
				statusMessageLabel.setText("Ready");
			}

		}

		class ResultPanel extends DiscMotifPanel {

			public ResultPanel() {
				super("Result");
				outputTextArea.setVisible(true);
				outputTextArea.setEditable(false);
				outputTextArea.setLineWrap(true);
				JScrollPane outputTextAreaScrollPane = new JScrollPane(outputTextArea);
				
				add(outputTextAreaScrollPane);

			}

		}

		class ReportPanel extends JPanel {
			public ReportPanel(JFileChooser chooser) {
				this.chooser = chooser;
				generateReportButton.setEnabled(false);
				generateReportButton.addActionListener(new GenerateReportAction(this.getParent()));
				sendToSTAMPButton.setEnabled(false);
				sendToSTAMPButton.addActionListener(new sendToSTAMPAction(this.getParent()));
				add(generateReportButton);
				add(sendToSTAMPButton);
			}


			class GenerateReportAction implements ActionListener {
				public GenerateReportAction(Container container) {
					this.parent = container;
				}
				@Override
				public void actionPerformed(ActionEvent event) {

					String reportDirectory = JOptionPane.showInputDialog(parent, "Please specify a directory to save the reports", chooser.getCurrentDirectory().toString() + File.separatorChar + "discmotif_reports");
					if (reportDirectory == null) {
						return;
					}
					File reportDir;

					try {
						reportDir = new File(reportDirectory);
					} catch (Exception e) {
						JOptionPane.showMessageDialog(parent, "Error with the report directory specified! Please set another one", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					}

					if (reportDir.exists()) {
						int selection = JOptionPane.showConfirmDialog(parent, "Warning: the directory you specified already exists. Overwrite previous report files?", "Overwrite?", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
						if (selection != JOptionPane.YES_OPTION) {
							return;
						}
					} else {
						try {			
							boolean mkdirSuccess = reportDir.mkdirs();
							if (!mkdirSuccess) {
								JOptionPane.showMessageDialog(parent, "Error: the directory you specified can not be created!", "Error", JOptionPane.ERROR_MESSAGE);
								return;
							}
						} catch (Exception e) {
							JOptionPane.showMessageDialog(parent, "Error while making directories! Please set another directory", "Error", JOptionPane.ERROR_MESSAGE);
							return;
						}
					}
					(new File(reportDirectory + File.separatorChar + "images")).mkdir();
					try {
						writeHtmlReport(reportDirectory);
					} catch (Exception e) {
						JOptionPane.showMessageDialog(parent, "Error while generating html report: " + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
						return;
					}

				}
				private Container parent;
			}

			class sendToSTAMPAction implements ActionListener {
				public sendToSTAMPAction(Container container) {
					parent = container;
				}
				public void actionPerformed(ActionEvent event) {
					BufferedReader rd;

					//int selection = JOptionPane.showInputDialog(parent, "<html>Clicking Yes below will send the discriminative motifs discovered by this program to STAMP, <br>which is an online tool for matching query motifs with known motifs in the TRANSFAC database. <br>Please allow several seconds for the result to be generated and loaded in your default browser. <br>Continue?</html>", "Submit to STAMP", JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE);
					String selectedDatabase = (String)JOptionPane.showInputDialog(parent, 
							"<html>Please select a database, and clicking \"OK\" below will send the discriminative <br>"+
							"motifs discovered by this program to STAMP, which is an online tool for matching <br>"+
							"query motifs with known motifs in the selected database. Please allow several seconds<br>"+
							"for the result to be generated and loaded in your default browser. <br>"+
							"Available databases:", "Submit to STAMP", JOptionPane.PLAIN_MESSAGE, null, possibleDatabases, "TRANSFAC v11.3");
					String queryDatabaseCode = "";
					if ((selectedDatabase == null) || (selectedDatabase.length()==0)) {
						return;
					} else {
						for (int i = 0; i < possibleDatabases.length; i++) {
							if (selectedDatabase.equals(possibleDatabases[i])) {
								queryDatabaseCode = possibleDatabaseCodes[i];
								break;
							}
						}
					}
					try {
						String toBeSent = URLEncoder.encode("input", "UTF-8") + "=" + URLEncoder.encode(outputStringBuilder.toString(), "UTF-8");
						toBeSent += "&" + URLEncoder.encode("match_db", "UTF-8") + "=" + URLEncoder.encode(queryDatabaseCode, "UTF-8");
						toBeSent += "&" + URLEncoder.encode("metric", "UTF-8") + "=" + URLEncoder.encode("PCC","UTF-8");
						toBeSent += "&" + URLEncoder.encode("align", "UTF-8") + "=" + URLEncoder.encode("SWU", "UTF-8");
						toBeSent += "&" + URLEncoder.encode("mulalign", "UTF-8") + "=" + URLEncoder.encode("IR", "UTF-8");
						toBeSent += "&" + URLEncoder.encode("tree", "UTF-8") + "=" + URLEncoder.encode("UPGMA", "UTF-8");
						toBeSent += "&" + URLEncoder.encode("num_hits", "UTF-8") + "=" + URLEncoder.encode("5", "UTF-8");
						toBeSent += "&" + URLEncoder.encode("mtrim", "UTF-8") + "=" + URLEncoder.encode("checked", "UTF-8");
						toBeSent += "&" + URLEncoder.encode("ICbox", "UTF-8") + "=" + URLEncoder.encode("0.4", "UTF-8");
						URL url = new URL(STAMP_URL);
						URLConnection conn = url.openConnection();
						conn.setDoOutput(true);
						OutputStreamWriter wr = new OutputStreamWriter(conn.getOutputStream());
						wr.write(toBeSent);
						wr.flush();

						rd = new BufferedReader(new InputStreamReader(conn.getInputStream()));
					} catch (Exception e) {
						JOptionPane.showMessageDialog(parent, "Failed to connect to STAMP server: " + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
						return;
					}

					Pattern p = Pattern.compile("http://.*html");
					String resultUrl = null;
					String line;
					try {
						while ((line = rd.readLine()) != null) {
							Matcher m = p.matcher(line);
							if (m.find()) {
								resultUrl = m.group(0);
								break;
							}
						}
					} catch (IOException e) {
						JOptionPane.showMessageDialog(parent, "Failed to run STAMP on server: " + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
						return;
					}
					if (resultUrl == null) {
						JOptionPane.showMessageDialog(parent, "Failed to parse the URL for STAMP result! ", "Error", JOptionPane.ERROR_MESSAGE);
						return;
					}
					String[] parts = resultUrl.split(":");
					try {
						java.awt.Desktop desktop = java.awt.Desktop.getDesktop();
						java.net.URI uri = new URI(parts[0], parts[1], null);  
						desktop.browse(uri);
					} catch (Exception e) {
						JOptionPane.showMessageDialog(parent, "Failed to open the STAMP result page in your default browser: " + e.toString(), "Error", JOptionPane.ERROR_MESSAGE);
						return;
					}

				}
				private Container parent;
			}


			private JFileChooser chooser;
		}

		class MainFrame extends JFrame {
			public MainFrame(JFileChooser chooser) {
				setTitle("DECOD: Fast and Accurate Discriminative Motif Finding");
				setSize(DEFAULT_WIDTH, DEFAULT_HEIGHT);
				setLocationByPlatform(true);

				GridBagLayout layout = new GridBagLayout();
				setLayout(layout);

				SequencePanel seqPanel = new SequencePanel(chooser);
				ParameterPanel paramPanel = new ParameterPanel(chooser);
				ActionPanel actionPanel = new ActionPanel(paramPanel.getOutputFilenameTextField(), this);
				ResultPanel resultPanel = new ResultPanel();
				StatusPanel statusPanel = new StatusPanel();
				ReportPanel reportPanel = new ReportPanel(chooser);

				//GUI for the sequence panel
				Font font = new Font("System", Font.BOLD, 14);
				JLabel titleLabel = new JLabel("DECOD: Deconvolved Discriminative Motif Finder " + VERSION_NUMBER);
				titleLabel.setFont(font);
				add(titleLabel, new GBCHelper(0,0,1,1).setAnchor(GBCHelper.WEST).setInsets(5,15,10,10));
				add(seqPanel, new GBCHelper(0,1,1,1).setAnchor(GBCHelper.WEST).setInsets(0,10,0,10));
				add(paramPanel, new GBCHelper(0,2,1,1).setWeight(100,0).setAnchor(GBCHelper.WEST).setInsets(0,10,0,10));
				add(actionPanel, new GBCHelper(0,3,1,1).setWeight(100,0));
				add(statusPanel, new GBCHelper(0,4,1,1).setWeight(100,0));
				add(resultPanel, new GBCHelper(0,5,1,1).setWeight(100,0).setInsets(0,5,0,5));
				add(reportPanel, new GBCHelper(0,6,1,1).setWeight(100,0));
			}
			private static final int DEFAULT_WIDTH = 600;
			private static final int DEFAULT_HEIGHT = 700;
		}
		//********BEGIN MAIN********//
		if (args.length>0 && args[0].equals("-nogui")) {
			commandLineMode = true;
			core(args);
		} 
		
		else if ((args.length > 0) && (!args[0].equals("-nogui")) && (!args[0].equals("-help"))) {
			commandLineMode = false;
			setupGlobalParameters(args);
			try {
				posSeqFilename = new File(posSeqFilename).getCanonicalPath();
			} catch (Exception e) {
				System.out.println("Error: positive sequence filename " + posSeqFilename + " is invalid!");
				e.printStackTrace();
				System.exit(1);
			}
			try {
				negSeqFilename = new File(negSeqFilename).getCanonicalPath();
			} catch (Exception e) {
				System.out.println("Error: negative sequence filename " + negSeqFilename + " is invalid!");
				e.printStackTrace();
				System.exit(1);
				
			}
			if (OUTPUT_FILENAME != null) {
				try {
					OUTPUT_FILENAME = new File(OUTPUT_FILENAME).getCanonicalPath();
				} catch (Exception e) {
					System.out.println("Error: output filename " + OUTPUT_FILENAME + " is invalid!");
					e.printStackTrace();
					System.exit(1);

				}
			}
			EventQueue.invokeLater(new Runnable()
			{
				public void run() {
					JFileChooser chooser;
					chooser = new JFileChooser();
					chooser.setCurrentDirectory(new File("."));
					MainFrame frame = new MainFrame(chooser);
					frame.setResizable(false);
					frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
					frame.setVisible(true);
				}

			});
		}
		
		else if (args.length == 0) {
			EventQueue.invokeLater(new Runnable()
			{
				public void run() {

					JFileChooser chooser;
					chooser = new JFileChooser();
					chooser.setCurrentDirectory(new File("."));
					MainFrame frame = new MainFrame(chooser);
					frame.setResizable(false);
					frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
					frame.setVisible(true);
				}

			});
		}


		else {
			commandLineMode = true;
			System.out.println("Error: '-nogui' must be explicitly specified for command-line mode.");
			String[] defaultArgs = {"-nogui", "-help"};
			core(defaultArgs);
		}
		
	
	}
}

class Comparer implements Comparator<Map.Entry<String, Double>> {
	public int compare(Map.Entry<String, Double> e1, Map.Entry<String, Double> e2) {
		Double i1 = (Double) e1.getValue();
		Double i2 = (Double) e2.getValue();
		return i2.compareTo(i1);
	}
}


class AbsComparer implements Comparator<Map.Entry<String, double[]>> {
	public int compare(Map.Entry<String, double[]> e1, Map.Entry<String, double[]> e2) {
		Double i1 = Math.abs((Double) e1.getValue()[1] - e1.getValue()[0]);
		Double i2 = Math.abs((Double) e2.getValue()[1] - e2.getValue()[0]);
		return i2.compareTo(i1);
	}
}

@SuppressWarnings("serial")
class GBCHelper extends GridBagConstraints
{

	public GBCHelper(int x, int y)
	{
		this.gridx = x;
		this.gridy = y;
	}

	public GBCHelper(int x, int y, int width, int height)
	{
		this.gridx = x;
		this.gridy = y;
		this.gridwidth = width;
		this.gridheight = height;
	}

	public GBCHelper setAnchor(int anchor)
	{
		this.anchor = anchor;
		return this;
	}


	public GBCHelper setWeight(double weightx, double weighty)
	{
		this.weightx = weightx;
		this.weighty = weighty;
		return this;
	}

	public GBCHelper setInsets(int top, int left, int bottom, int right)
	{
		this.insets = new Insets(top, left, bottom, right);
		return this;
	}

}

