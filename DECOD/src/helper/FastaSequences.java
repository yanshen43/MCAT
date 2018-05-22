package helper;

import java.util.ArrayList;

public class FastaSequences {
	private ArrayList<String> ids;
	private ArrayList<String> seqs;
	public FastaSequences(ArrayList<String> ids, ArrayList<String> seqs) {
		this.ids = ids;
		this.seqs = seqs;
	}
	public ArrayList<String> getIds() {
		return this.ids;
	}
	public ArrayList<String> getSeqs() {
		return this.seqs;
	}

}
