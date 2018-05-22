package helper;

import java.util.ArrayList;
public class Result {
	private double score;
	private ArrayList<String> sites;
	private Pwm theta;
	private double[] pEsts;
	private Double pvalue;
	public Result() {
		this.score = 0;
		this.sites = new ArrayList<String>();
		this.pEsts = new double[2];
	}
	public Result(Result a) {
		this.score = a.getScore();
		this.sites = a.getSites();
		this.theta = a.getPwm();
		this.pEsts = a.getPEsts();
	}
	public Result(double score, ArrayList<String> sites, double[] pEsts) {
		this.score = score;
		this.sites = sites;
		this.pEsts = pEsts;
	}
	public Result(double score, ArrayList<String> sites, Pwm theta, double[] pEsts) {
		this.score = score;
		this.sites = sites;
		this.theta = theta;
		this.pEsts = pEsts;
	}
	public Result(double score, Pwm theta, double[] pEsts) {
		this.score = score;
		this.sites = null;
		this.theta = theta;
		this.pEsts = pEsts;
	}
	public void setPValue(Double pvalue) {
		this.pvalue = pvalue;
	}
	
	public Double getPValue() {
		return this.pvalue;
	}
	public double getScore() {
		return this.score;
	}
	public ArrayList<String> getSites() {
		return new ArrayList<String>(this.sites);
	}
	public Pwm getPwm() {
		return new Pwm(this.theta);
	}
	public double getPEstPos() {
		return pEsts[1];
	}
	public double getPEstNeg() {
		return pEsts[0];
	}
	public double[] getPEsts() {
		return pEsts;
	}
	public void setPEsts(double[] pEsts) {
		this.pEsts = pEsts;
	}
	public String toString() {
		StringBuilder sb = new StringBuilder();
		sb.append("score="+score+"\n");
		sb.append("theta="+theta+"\n");
		sb.append("pEstPos="+pEsts[1]+" pEstNeg= "+ pEsts[0] + "\n");
		return sb.toString();
	}
}
