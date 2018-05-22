package helper;

public class MotifInstance {
	int position;
	double score;
	String kMer;
	String id;
	public MotifInstance(String id, String kMer, int pos, double score) {
		this.id = id;
		this.kMer = kMer;		
		this.position = pos;
		this.score = score;
	}
	public int getPosition() {
		return position;
	}
	public String getId() {
		return this.id;
	}
	public double getScore() {
		return score;
	}
	public String getKMer() {
		return kMer;
	}
	public String toString() {
		return ">" + id + "\t" + position + "\t" + kMer + "\t" + String.format("%.4f",score);
	}	
}

