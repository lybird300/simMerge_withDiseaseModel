package seqSIM;

import java.io.Serializable;

public class Subject implements Comparable<Subject>, Serializable{
	private static final long serialVersionUID = 1L;
	private int ID;
	private int chromA;
	private int chromB;
	private int carrierOrNot;//0-non-carrier, 1-carrier
	public Subject(int id, Integer chrom1, Integer chrom2, int carryStatus) {
		this.ID = id;
		this.chromA = chrom1.intValue();
		this.chromB = chrom2.intValue();
		carrierOrNot = carryStatus;
	}
	public int isCarrier() {
		return carrierOrNot;
	}
	public void setCarrier(int carried) {
		this.carrierOrNot = carried;
	}
	public int getID() {
		return ID;
	}
	public void setID(int id) {
		ID = id;
	}
	public int getChromA() {
		return chromA;
	}
	public void setChromA(int chromA) {
		this.chromA = chromA;
	}
	public int getChromB() {
		return chromB;
	}
	public void setChromB(int chromB) {
		this.chromB = chromB;
	}
	public int compareTo(Subject theOtherSubject) {
		return this.getID() - theOtherSubject.getID();
	}
}

