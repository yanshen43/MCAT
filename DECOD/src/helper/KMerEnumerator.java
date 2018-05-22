package helper;
import java.io.FileOutputStream;
import java.io.PrintWriter;


public class KMerEnumerator {
	private static String kMer;
	private static int width;

	public KMerEnumerator(int w){
		kMer = "";
		width = w;
	}
	
	public void reset() {
		kMer = "";
	}
	
	private boolean tryAdd(int pos) {
		char c = kMer.charAt(pos);
		if (c < '3') {
			kMer = kMer.substring(0, pos) + (++c) + kMer.substring(pos+1);
			return false;
		} else if (pos==0) {
			return true;
		} else {
			kMer = kMer.substring(0, pos) + '0' + kMer.substring(pos+1);
			return tryAdd(pos-1);
		}
	}
	
	
	public String getNext() {
		if (kMer == null) {
			return null;
		}
		
		if (kMer == "") {
			StringBuilder sb = new StringBuilder();
			for (int i = 0; i < width; i++) {
				sb.append("0");
			}
			kMer = sb.toString();
			return kMer;
		}
		
		boolean end = tryAdd(width-1);
		
		if (end) {
			return null;
		}
		
		return kMer;
		
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		KMerEnumerator s = new KMerEnumerator(8);
		
		String kmer;
		PrintWriter pw = null;
		try {
		pw = new PrintWriter(new FileOutputStream("test.txt"));
		} catch (Exception e) {
			
		}
		while ((kmer = s.getNext()) != null) {
			pw.println(kmer);
			
		}
		pw.close();
		
	}


}
