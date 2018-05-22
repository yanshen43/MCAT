package helper;
import java.awt.Color;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics2D;
import java.awt.font.GlyphVector;
import java.awt.geom.AffineTransform;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
import java.util.Map;
import java.util.Scanner;
import java.util.Map.Entry;

import javax.imageio.ImageIO;

/**
 * 
 */

/**
 * @author Shan
 *
 */
public class Pwm {
	private int width;
	private double[][] matrix;
	private ArrayList<String> siteSets;
	public double[][] scoringMatrix;

	/**
	 * Create a Pwm object from an input file including probabilities of each site and base.
	 * No pseudocount added this way
	 * @param filename
	 */
	public Pwm(String filename) {
		String line;
		Scanner inputStream = null;
		ArrayList<ArrayList<Double>> matrixList = new ArrayList<ArrayList<Double>>();
		try {
			inputStream=new Scanner(new FileInputStream(filename));  
		} catch (FileNotFoundException e) {
			System.out.println("Fatal error: file "+filename+" was not found or could not be opened!");
			System.exit(1);
		}
		while (inputStream.hasNextLine()) {
			line = inputStream.nextLine();
			line.trim();
			if (line.startsWith("#")) continue;
			if (line.startsWith(">")) continue;
			if (line.startsWith("Site set")) break;
			if (line.startsWith("A") || line.startsWith("T") || line.startsWith("C") || line.startsWith("G")) {
				line = line.substring(3, line.length()-1);
			}
			
			String[] tokens = line.split(" ");
			
			ArrayList<Double> row = new ArrayList<Double>();
			for (String element:tokens) {
				row.add(Double.parseDouble(element));
			}
			matrixList.add(row);
		}//while
		if (matrixList.size()!=4) {
			System.out.println("Fatal error: PWM matrix file " + filename +"contains less or more than 4 lines!");
			System.exit(1);
		}
		this.width = matrixList.get(0).size();
		
		matrix = new double[4][width];
		for (int i=0; i<4; i++) {
			for (int j=0; j<width; j++) {
				matrix[i][j] = matrixList.get(i).get(j);
			}
		}
		scoringMatrix = new double[4][width];
	}
	
	
	/**
	 * Create a Pwm object from a list of sequences. Pseudocount of 0.1 is added
	 * The sequences must be of equal length
	 * @param siteSets
	 */
	public Pwm(ArrayList<String> siteSets){
		this.siteSets = new ArrayList<String>(siteSets);
		String oneS = siteSets.get(0);
		this.width = oneS.length();
		this.matrix = new double[4][width];
		for (String s: siteSets) {
			if (s.length()!=this.width) {
				System.out.println("Fatal error: the site-sets thrown in for generating PWM have illegal lengths!");
				System.out.println("current kmer  = "+s);
				System.out.println("All kmers = " + siteSets);
				System.out.println("Specified width="+this.width);
				System.exit(1);
			}
			char[] letters = s.toCharArray();
			for (int j = 0; j < width; j++) {
				matrix[letters[j]-48][j]++; //"0"-"9" corresponds to ascii 48-57
			}
		}
		for (int i=0; i<width; i++) {
			for (int j=0; j<4; j++) {
				matrix[j][i] = (matrix[j][i]+0.1)/(siteSets.size()+0.4); //add pseudocounts
			}
		}
		scoringMatrix = new double[4][width];
	}
	public Pwm(double[][] matrix) {
		this.matrix = matrix.clone();
		this.width = matrix[0].length;
		this.scoringMatrix = new double[4][width];
	}

	public Pwm(Pwm anotherPwm) {
		this.width = anotherPwm.getWidth();
		this.matrix = anotherPwm.getMatrix();
		this.siteSets = anotherPwm.getSiteSets();
		this.scoringMatrix = anotherPwm.getScoringMatrix();
	}

	
	public double[][] getMatrixWithoutPseudoCounts() {
		if (this.siteSets == null) {
			return this.getMatrix();
		}
		double[][] m = new double[4][this.getWidth()];
		for (String s: siteSets) {
			char[] letters = s.toCharArray();
			for (int j = 0; j < width; j++) {
				m[letters[j]-48][j]++; //"0"-"9" corresponds to ascii 48-57
			}
		}
		for (int i=0; i<width; i++) {
			for (int j=0; j<4; j++) {
				m[j][i] = (m[j][i])/(siteSets.size());
			}
		}
		return m;
	}
	
	
	public ArrayList<String> getSiteSets() {
		return new ArrayList<String>(this.siteSets);
	}
	public void setMatrix(double[][] m) {
		this.matrix = m;
	}
	public Integer getWidth(){
		return width;
	}
	public double getElement(Integer m, Integer n){
		return matrix[m][n];
	}
	
	public double[][] getMatrix() {
		double[][] copyMatrix = new double[4][width];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < width; j++) {
				copyMatrix[i][j] = this.matrix[i][j];
			}
		}
		
		return copyMatrix;
	}
	
	public double[][] getScoringMatrix() {
		if (this.scoringMatrix == null) {
			return null;
		}
		double[][] copyMatrix = new double[4][width];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < width; j++) {
				copyMatrix[i][j] = this.scoringMatrix[i][j];
			}
		}
		
		return copyMatrix;
	}
	
	public void convertToScoringMatrix(Pwm background) {
		for (int i = 0; i < 4; i++) {
			for (int j=0; j < width; j++) {
				this.scoringMatrix[i][j] = Math.log( this.matrix[i][j] / background.getElement(i,0));
			}
		}
		
	}
	public String toString() {
		StringBuilder sb = new StringBuilder();
		double[][] m = this.getMatrixWithoutPseudoCounts();
		for (int i=0; i<4; i++) {
			switch (i) {
				case 0: sb.append("A [");break;
				case 1: sb.append("C [");break;
				case 2: sb.append("G [");break;
				case 3: sb.append("T [");break;
			}
			for (int j=0; j<width; j++) {
				if (j!=width-1) {
					sb.append(String.format("%.4f ",Math.round(m[i][j]*10000)/10000.0));
					//sb.append(m[i][j]+" ");
				} else {
					sb.append(String.format("%.4f",Math.round(m[i][j]*10000)/10000.0));
					//sb.append(m[i][j]);
				}
			}
			if (i!=3) {sb.append("]\n");
			} else {
				sb.append("]");
			}
		}
		return sb.toString();
	}


	public String toHtml() {
		StringBuilder sb = new StringBuilder();
		sb.append("<font face=\"Courier, Courier New\">");
		double[][] m = this.getMatrixWithoutPseudoCounts();
		for (int i=0; i<4; i++) {
			switch (i) {
				case 0: sb.append("A [");break;
				case 1: sb.append("C [");break;
				case 2: sb.append("G [");break;
				case 3: sb.append("T [");break;
			}
			for (int j=0; j<width; j++) {
				if (j != width-1) {
					sb.append(String.format("%.2f ", m[i][j]));
				} else {
					sb.append(String.format("%.2f", m[i][j]));
				}
			}
			if (i!=3) {sb.append("]<br>");
			} else {
				sb.append("]");
			}
		}
		sb.append("</font>");
		return sb.toString();
	}
	
	public double[] getColumn(int col) {
		double[] column = new double[4];
		for (int i=0; i<4; i++) {
			column[i] = matrix[i][col];
		}
		return column;
	}

	public double emit(String s) {
		if ((this.width!=1) && (s.length()!=this.width)) {
			System.out.println("Fatal error: length of the sequence "+s+" does not match the width of the PWM!");
			System.exit(1);
		}
		double prob = 1;
		char[] letters = s.toCharArray();
		for (int i=0; i<s.length(); i++) {
			if (this.width == 1) {
				prob *= this.matrix[letters[i]-48][0]; //Background model
			}
			else {
				prob *=  this.matrix[(letters[i]-48)][i];
			}
		} 
		return prob;
	}
	
	public double score(String s) {
		if ((this.width!=1) && (s.length()!=this.width)) {
			System.out.println("Fatal error: length of the sequence "+s+" does not match the width of the PWM!");
			System.exit(1);
		}
		double score = 0;
		char[] letters = s.toCharArray();
		for (int i=0; i<s.length(); i++) {
			score += this.scoringMatrix[(letters[i]-48)][i];
		}
		return score;
	}
	
	public double partEmit(int[] positions, String s) { //positions: one-based!!
		if (positions.length != s.length()) {
			System.out.println("Fatal error: length of string "+s+" does not match the length of the given positions!");
			System.exit(1);
		}
		char[] letters = s.toCharArray();
		double prob = 1;
		for (int i = 0; i < s.length(); i++) {
			//System.out.println(i+" "+positions[i]+" "+(letters[i]-48));
			prob *= this.matrix[letters[i]-48][positions[i]-1];
		}
		return prob;
	}
	
	
	public double partEmit(int start, int end, String s) { //one-based!!
		if ((start<1) || (end < start) || (end > this.width) || (end-start+1 != s.length())) {
			System.out.println("Fatal error: length of the sequence "+s+" does not match the width of the PWM starting at "+start+" and ending at "+end+"!");
			System.exit(1);
		}
		double prob = 1;
		char[] letters = s.toCharArray();
		for (int i = 0; i<s.length(); i++) {
			prob*=this.matrix[letters[i]-48][i+start-1];
		}
		return prob;
	}
	
		
	public String generateSequence() {
		StringBuffer sb = new StringBuffer();
		SecureRandom generator = new SecureRandom();
		for (int i=0; i<width; i++) {
			double[] column = this.getColumn(i);
			double[] cdf = new double[4];
			for (int j=0; j<4; j++) {
				for (int l=0; l<=j; l++) {
					cdf[j]+=column[l];
				}
			}
			
			double r = generator.nextDouble();
			for (int j=0;j<4; j++) {
				if (r<=cdf[j]) {
					sb.append(j);
					break;
				}
			}
		}
		return sb.toString();
	}
	
	public double[][] minus(Pwm another) {
		if (this.getWidth()!=another.getWidth()) {
			System.out.println("Fatal error: The two PWMs have different lengths!");
			System.exit(1);
		}
		double[][] diff = new double[4][this.width];
		for (int i=0; i<4; i++) {
			for (int j=0; j<this.width; j++) {
				diff[i][j] = this.getElement(i,j)-another.getElement(i,j);
			}
		}
		return diff;
	}
	
	
	
	public double KLDivergence(Pwm anotherPwm) {
		double s = 0;
		if (anotherPwm.getWidth() != this.getWidth()) {
			System.out.println("Fatal error: The widths of the two motifs being calculated for KL divergence are different!");
			System.exit(1);
		}
		double[][] m1 = this.getMatrixWithoutPseudoCounts();
		double[][] m2 = anotherPwm.getMatrixWithoutPseudoCounts();
		for (int j = 0; j < this.getWidth(); j++) {
			for (int i=0; i<4; i++) {
				if (m2[i][j]==0 && m1[i][j]!=0) {
					m2[i][j]=0.005;
				}
				if (m1[i][j]==0 && m2[i][j]!=0) {
					m1[i][j]=0.005;
				}
				if (m1[i][j]!=m2[i][j]) {
					s += (m1[i][j]-m2[i][j])*Math.log(m1[i][j]/m2[i][j]);
				}
			}
		}
		return s;
		
	}
	
	public void convertToLogo(String filename) {
		int height = 200;
		int singleLetterWidth = 60;
		
		double[][] m = this.getMatrixWithoutPseudoCounts();
		
		double[] columnICs = new double[this.getWidth()];
		int[][] letterHeights = new int[4][this.getWidth()];
		
		for (int i = 0; i < this.getWidth(); i++) {
			for (int j = 0; j < 4; j++) {
				if (m[j][i] != 0) {
					columnICs[i] -= m[j][i] * Math.log(m[j][i])/Math.log(2);
				}
			}
			columnICs[i] = 2 - columnICs[i];
		}
		
		
		
		for (int i = 0; i < this.getWidth(); i++) {
			for (int j = 0; j < 4; j++) {
				if (m[j][i] == 0) {
					letterHeights[j][i] = 0;
				} else {
					letterHeights[j][i] = (int) Math.round(m[j][i] * columnICs[i]/2 * height);
				}
			}
		}

		ArrayList<List<Map.Entry<Integer, Integer>>> sortedLetterHeights = new ArrayList<List<Map.Entry<Integer, Integer>>>(); // cumulative heights!
		ArrayList<List<Map.Entry<Integer, Integer>>> sortedCumulativeLetterHeights = new ArrayList<List<Map.Entry<Integer, Integer>>>(); // cumulative heights!
		ArrayList<List<Map.Entry<Integer, Integer>>> sortedCumulativeLetterBaselines = new ArrayList<List<Map.Entry<Integer, Integer>>>(); // cumulative heights!
		for (int i = 0; i < this.getWidth(); i++) {
			Hashtable<Integer, Integer> map = new Hashtable<Integer, Integer>();
			for (int j = 0; j < 4; j++) {
				map.put(j, letterHeights[j][i]);
			}
			Hashtable<Integer, Integer> map2 = new Hashtable<Integer, Integer>(map);
			Hashtable<Integer, Integer> map3 = new Hashtable<Integer, Integer>(map);
			
			List<Map.Entry<Integer, Integer>> list = new ArrayList<Map.Entry<Integer, Integer>>(map.entrySet());
			List<Map.Entry<Integer, Integer>> list2 = new ArrayList<Map.Entry<Integer, Integer>>(map2.entrySet());
			List<Map.Entry<Integer, Integer>> list3 = new ArrayList<Map.Entry<Integer, Integer>>(map3.entrySet());
			
			Comparer c = new Comparer();
			Collections.sort(list, c);
			Collections.sort(list2, c);
			Collections.sort(list3, c);
			
			sortedLetterHeights.add(list);
			sortedCumulativeLetterHeights.add(list2);
			sortedCumulativeLetterBaselines.add(list3);
		
		}

		/*System.out.println("Letter heights");
		for (List<Entry<Integer, Integer>> col : sortedLetterHeights) {
			for (Entry<Integer, Integer> letter : col) {
				int i = letter.getKey();
				int h = letter.getValue();
				System.out.println(i+" "+h);
			}
			System.out.println("***");
		}		
		*/
		//Calculate cumulative letter heights up to the current letter
		//ArrayList<List<Map.Entry<Integer, Integer>>> sortedLetterCumulativeHeights = new ArrayList<List<Map.Entry<Integer, Integer>>>(sortedLetterHeights);
		for (List<Entry<Integer, Integer>> col : sortedCumulativeLetterHeights) {
			int cumulative = 0;
			for (Entry<Integer, Integer> letter : col) {
				int h = letter.getValue();
				cumulative += h;
				letter.setValue(cumulative);
			}
		}

		
		//ArrayList<List<Map.Entry<Integer, Integer>>> sortedLetterCumulativeBaselines = new ArrayList<List<Map.Entry<Integer, Integer>>>(sortedLetterHeights);
		for (List<Entry<Integer, Integer>> col : sortedCumulativeLetterBaselines) {
			int cumulative = 0;
			for (Entry<Integer, Integer> letter : col) {
				int h = letter.getValue();
				letter.setValue(cumulative);
				cumulative += h;
			}
		}
		
		/*System.out.println("Letter heights (cumulative)");
		for (List<Entry<Integer, Integer>> col : sortedCumulativeLetterHeights) {
			for (Entry<Integer, Integer> letter : col) {
				int i = letter.getKey();
				int h = letter.getValue();
				System.out.println(i+" "+h);
			}
			System.out.println("***");
		}*/
		
		/*System.out.println("Letter baselines");
		for (List<Entry<Integer, Integer>> col : sortedCumulativeLetterBaselines) {
			for (Entry<Integer, Integer> letter : col) {
				int i = letter.getKey();
				int h = letter.getValue();
				System.out.println(i+" "+h);
			}
			System.out.println("***");
		}*/
		
		int colNumber = 0;
		
		BufferedImage bufferedImage = new BufferedImage((singleLetterWidth+5) * this.getWidth(), height, BufferedImage.TYPE_INT_RGB);
		Graphics2D g2d = bufferedImage.createGraphics();
		g2d.setColor(Color.white);
		g2d.fillRect(0,0,(singleLetterWidth+5) * this.getWidth(),height);
		
		for (List<Entry<Integer, Integer>> col : sortedLetterHeights) {
			int rowNumber = 0;
			for (Entry<Integer, Integer> letter : col) {
				int i = letter.getKey();
				int h = letter.getValue();
				String s = "";
				Color c = Color.green;
				
				switch (i) {
					case 0: {s = "A"; c = Color.green; break;}
					case 1: {s = "C"; c = Color.blue; break;}
					case 2: {s = "G"; c = Color.yellow; break;}
					case 3: {s = "T"; c = Color.red; break;}
				}
			
				Font font = new Font("SansSerif", Font.BOLD, 150);
				g2d.setFont(font);
				//Rectangle2D fontDim = font.getStringBounds(s, g2d.getFontRenderContext());
				FontMetrics fmOld = g2d.getFontMetrics(font);
				int letterWidth = fmOld.stringWidth(s);
				//int letterHeight = 88;//corresponding to font point size 120
				GlyphVector gv = font.createGlyphVector(g2d.getFontRenderContext(),s);
				Rectangle2D bounds = gv.getGlyphMetrics(0).getBounds2D();
				int letterHeight = (int)bounds.getHeight();
				
				
				font = font.deriveFont(AffineTransform.getScaleInstance((double)singleLetterWidth/letterWidth, (double)h/letterHeight)); 
				g2d.setFont(font);
				g2d.setPaint(c);
				int thisBaseline = sortedCumulativeLetterBaselines.get(colNumber).get(rowNumber).getValue();
				//System.out.println(s+" "+(colNumber*singleLetterWidth)+" "+(height-thisBaseline)+" "+g2d.getFontMetrics().getAscent());
				g2d.drawString(s, colNumber*singleLetterWidth, height - thisBaseline);
				rowNumber++;
				
			}
			colNumber++;
		}
		try {
			ImageIO.write(bufferedImage, "PNG", new File(filename));
		} catch (Exception e) {
			e.printStackTrace();
		}
		

	}
	
	
	public Pwm getReverseComplement() {
		double[][] matrix = this.getMatrix();
		double[][] reverseMatrix = new double[matrix.length][matrix[0].length];
		for (int i=0; i<matrix[0].length; i++) {
			for(int l=0; l<matrix.length; l++) {
				reverseMatrix[matrix.length-l-1][matrix[0].length-i-1] = matrix[l][i];
			}
		}
		return new Pwm(reverseMatrix);
	}
	
	public static void main(String args[]) {
		
		Pwm theta1 = new Pwm("./data/simulation/width6/100seqs/one_motif/1.15/theta.pwm");
		Pwm theta2 = new Pwm("./data/simulation/width6/100seqs/one_motif/0.64/theta.pwm");
		
		System.out.println(theta1.KLDivergence(theta2));
		

	}
	
	class Comparer implements Comparator<Map.Entry<Integer, Integer>> {
		public int compare(Map.Entry<Integer, Integer> e1, Map.Entry<Integer, Integer> e2) {
			Integer i1 = (Integer) e1.getValue();
			Integer i2 = (Integer) e2.getValue();
			return i1.compareTo(i2);
		}
	}

}
	