/*
 * Originally written as part of a larger project by Aarathi Raghuraman (raarathi@vt.edu)
 * Adapted as a stand alone tool by Jeff Robertson (thejar@vt.edu)
 */
import java.io.*;
import java.util.*;

public class IntuitionisticScoring
{
    private class Motif {
        public double[][] pwm;
        public String consensus;
        public int width;
        public void readPWM(String pwmFilename)
            throws FileNotFoundException
        {
            pwm = new double[width][4];
            Scanner in = new Scanner(new File(pwmFilename));
            for (int w = 0; w < width; w++) {
                for (int letter = 0; letter < 4; letter++) {
                    pwm[w][letter] = in.nextDouble();
                }
            }
        }
    }
	public double alpha = 0.01;
	private Motif mtf;
	private double[][] pwm;
	private String motif;
	private double[][] mem = new double[16][16];
	private double[][] nonMem = new double[16][16];
	private double[][] intuit = new double[16][16];
	private double[][] normIntuit = new double[16][16];
	private ArrayList<String> sequence;

    public static void main(String[] args)
        throws FileNotFoundException
    {
        String pwmFilename = args[0];
        String fastaFilename = args[1];
        int width = Integer.parseInt(args[2]);
        IntuitionisticScoring is = new IntuitionisticScoring(pwmFilename, width, fastaFilename);
        System.out.println(is.score());
    }

	public IntuitionisticScoring()
	{
		mtf = new Motif();
	}

	public IntuitionisticScoring(String pwmFilename, int width, String fileName)
        throws FileNotFoundException
	{
		mtf = new Motif();
        Motif m = mtf;
        mtf.width = width;
        mtf.readPWM(pwmFilename);
		pwm = m.pwm;
		motif = m.consensus;
		if(motif == null)
		{
			motif = "";
			for(int w=0; w<m.width; w++)
			{
				double  max = 0;
				int maxLetter = -1;
				for(int letter = 0; letter<4;letter++)
				{
					if(pwm[w][letter]>max)
					{
						max = pwm[w][letter];
						maxLetter = letter;
					}
				}
				motif += getIntToChar(maxLetter);
			}
		}

		for(int letter=0; letter<4; letter++)
		{
			for(int loc=0; loc<motif.length();loc++)
			{
				if(pwm[loc][letter] == 0)
				{
					pwm[loc][letter] = alpha;
				}
			}
		}
		m.consensus = motif;
		ArrayList<String> sequences = new ArrayList<String>();
        Scanner in = new Scanner(new File(fileName));
    
        String line = "";
        char c;
        String str = "";
        while(in.hasNext())
        {
            line = in.nextLine();
            if (line.length() > 0 && line.charAt(0) != '>' && line.charAt(0) != ' ')
            {
                    for (int i = 0; i < line.length(); i++)
                    {
                                c = line.charAt(i);
                                    str += c;
                    }
            }
            else
            {
                    if((str.trim()).length()>0)
                    {
                            sequences.add(str);
                            str = "";
                    }
            }
        }
		sequence = sequences;
	}

	public double score()
	{
		int n = motif.length();
		double scIntuit = 0.0;
		double maxSCIntuit = 0.0;
		double avg  = 0.0;
		for(int k =0; k<sequence.size(); k++)
		{
			int count = 0;
			for(int l =0; l<((sequence.get(k)).length())-n;l++)
			{
				scIntuit = 0.0;
				for(int i = l; i<l+n-1; i++)
				{
					for(int j = i+1; j<l+n; j++)
					{
						pos1pos2Score(i-l,j-l);
						scIntuit += normIntuit[getCharToInt(sequence.get(k).charAt(i))][
                            getCharToInt(sequence.get(k).charAt(j))];
					}
				}
				scIntuit /= ((n-1)*n)/2;
				maxSCIntuit = scIntuit>maxSCIntuit?scIntuit:maxSCIntuit;
			}
			avg += maxSCIntuit;
		}
		avg /= sequence.size();
		return avg;
	}

	public void pos1pos2Score(int pos1, int pos2)
	{
		double maxNonMem = 0.0;
		// Finds the value of membership, non-membership and the max 
        // non-membership for all possible (16) nucleotide combinations
		for(int letter = 0; letter<4; letter++)
		{
			for(int letter2 = 0; letter2<4; letter2++)
			{
				mem[letter][letter2] = membership(letter, letter2, pos1, pos2);
				nonMem[letter][letter2] = nonmembership(letter, letter2, pos1, pos2);
				if(maxNonMem < nonMem[letter][letter2])
				{
					maxNonMem = nonMem[letter][letter2];
				}
			}
		}
		
		double maxIntuit = 0.0;
		double minIntuit = 1.0;
		// Finding the score of every combination of nucleotides at position i, j
		for(int letter = 0; letter<4; letter++)
		{
			for(int letter2 = 0; letter2<4; letter2++)
			{
				intuit[letter][letter2] = mem[letter][letter2]*(maxNonMem-nonMem[letter][letter2]);
				if(maxIntuit < intuit[letter][letter2])
				{
					maxIntuit = intuit[letter][letter2];
				}
				if(minIntuit > intuit[letter][letter2])
				{
					minIntuit = intuit[letter][letter2];
				}
			}
		}

		// Calculates the normalized score for pos1, pos2 for every possible combination of nucleotides
		for(int letter = 0; letter<4; letter++)
		{
			for(int letter2 = 0; letter2<4; letter2++)
			{
				normIntuit[letter][letter2] = (intuit[letter][letter2]
                        - minIntuit)/(maxIntuit-minIntuit);
			}
		}
	}

	public double membership(int b1, int b2, int pos1, int pos2)
	{
		double pairProb = pwm[pos1][b1]*pwm[pos2][b2];
		double membershipProb = pairProb+((1-pairProb)*((pwm[pos1][b1]+pwm[pos2][b2])/2));
		return membershipProb;
	}

	public double nonmembership(int b1, int b2, int pos1, int pos2)
	{
		double membershipProb = membership(b1, b2, pos1, pos2);
		double nonmembershipProb = ((ICbp(b1, pos1)+ICbp(b2, pos2))/2)*(1-membershipProb);
		return nonmembershipProb;
	}

	public double ICbp(int b, int pos)
	{
		double val = (2+(pwm[pos][b]*(Math.log(pwm[pos][b])/Math.log(2))))/2;
		return val;
	}

	public int getCharToInt(char nucleotide)
	{
		switch(nucleotide)
		{
			case 'A':
				return 0;
			case 'a':
				return 0;
			case 'C':
				return 1;
			case 'c':
				return 1;
			case 'G':
				return 2;
			case 'g':
				return 2;
			case 'T':
				return 3;
			case 't':
				return 3;
			default:
				return -1;
		}
	}

	public char getIntToChar(int nucleotide)
	{
		switch(nucleotide)
		{
			case 0:
				return 'A';
			case 1:
				return 'C';
			case 2:
				return 'G';
			case 3:
				return 'T';
			default:
				return '*';
		}
	}
}
