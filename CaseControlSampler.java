package seqSIM;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;

public class CaseControlSampler /*implements Runnable*/ {
	private int grrIndex;
	private int bdrIndex;
	private int totalCase=0;
	private int totalCtrl=0;
	private String location;
	Set<Integer> usedChroms=null;

	public CaseControlSampler(File parentDir, int grrIndex, int bdrIndex, int repIndex) {
		this.grrIndex = grrIndex;
		this.bdrIndex = bdrIndex;
		if(ParamParser.createSingleSet){
			totalCase += ParamParser.num_ca[0];
			totalCtrl += ParamParser.num_co[0];
		}
		else{
			for(int subSetID = 0; subSetID < ParamParser.num_subsets; subSetID++){
				totalCase += ParamParser.num_ca[subSetID];
				totalCtrl += ParamParser.num_co[subSetID];
			}
		}
		/*At the 3rd level directory, create a folder for the current case-control panel and name the folder as:
		 * Case<numOfCase>_Control<numOfControl>_grr<valueOfGRR>*/
		String workDir = "SubsetCt"+ParamParser.num_subsets+"_GRR"+ParamParser.grr[grrIndex]+"_BDR"+ParamParser.bdr[bdrIndex]+File.separator+"Rep"+repIndex;
		File dir_3rd = new File(parentDir, workDir);
		if(!dir_3rd.exists()) dir_3rd.mkdirs();
		location = dir_3rd.getAbsolutePath();
		usedChroms = new HashSet<Integer>();
	}
	
	public String getLocation(){
		return location;
	}

	public List<DataSet> prepPanel(double cvMaf, List<Integer> chromsWithMut, List<Integer> chromsWithoutMut) { 		
	  	List<DataSet> allsets = new ArrayList<DataSet>();
	  	DataSet uniSet = assembCaseControlPanel(chromsWithMut, chromsWithoutMut, cvMaf);
	  	allsets.add(uniSet);
	  	if(ParamParser.createSingleSet){
	  		Util.outputFamFile(new File(location, ParamParser.outPrefix + "Set1.fam"), uniSet.caseList, uniSet.ctrlList);
	  	}
	  	else{
	  		Util.outputFamFile(new File(location, ParamParser.outPrefix + "Set" + uniSet.dataSetID + ".fam"), uniSet.caseList, uniSet.ctrlList);
		  	allsets.addAll(dividePanel(uniSet.caseList, uniSet.ctrlList));
	  	}
	  	return allsets;	  	
	}
	
	public int[][] getsubHapID(){
		int[][] subHapID = new int[totalCase + totalCtrl][2];
		File subInfoFile = null;
		if(ParamParser.createSingleSet)
			subInfoFile = new File(location, ParamParser.outPrefix + "singleSetSubInfo.csv");
		else
			subInfoFile = new File(location, ParamParser.outPrefix + "allSubInfo.csv");
		FileHandler readSubInfo = new FileHandler(subInfoFile);
		readSubInfo.openForRead();
		readSubInfo.nextLine();//header
		String nextLine = null;
		while((nextLine = readSubInfo.nextLine()) != null){
			String[] s = nextLine.split(",");
			int cellID = Integer.valueOf(s[0]) - 1;
			subHapID[cellID][0] = Integer.valueOf(s[1]);
			subHapID[cellID][1] = Integer.valueOf(s[2]);
		}
		readSubInfo.closeReader();
		return subHapID;
	}

	private List<DataSet> dividePanel(List<Subject> caseList, List<Subject> controlList) {
		if(caseList==null || controlList==null){
			System.out.println("A list parameter is null.");
			System.exit(0);
		}
		if(caseList.size()<=0 || controlList.size() <= 0){
			System.out.println("A list is empty.");
			System.exit(0);
		}
		if(ParamParser.num_subsets <= 0){
			System.out.println("The number of divided sets must be > 0");
			System.exit(0);
		}
		if(ParamParser.num_subsets == 1)
			return null;
		List<DataSet> subsets = new ArrayList<DataSet>();
		int caseStart = 0;
		int ctrlStart = 0;
		for(int subsetID = 0; subsetID < ParamParser.num_subsets; subsetID++){
			DataSet ds = new DataSet(subsetID+1, ParamParser.num_ca[subsetID], ParamParser.num_co[subsetID]);
			int caseEnd = caseStart + ParamParser.num_ca[subsetID];
			int ctrlEnd = ctrlStart + ParamParser.num_co[subsetID];
			ds.caseList.addAll(caseList.subList(caseStart, caseEnd));
			ds.ctrlList.addAll(controlList.subList(ctrlStart, ctrlEnd));
			caseStart = caseEnd;
			ctrlStart = ctrlEnd;
			Util.outputFamFile(new File(location, ParamParser.outPrefix + "Set" + ds.dataSetID + ".fam"), ds.caseList, ds.ctrlList);
			subsets.add(ds);
		}
		return subsets;
	}
	
	public int getNextIndex(List<Integer> listToCheck, int curIndex){
		if(curIndex == listToCheck.size()){
			Collections.shuffle(listToCheck);
			return 0;
		}
		return ++curIndex;
	}
	
	/**
	 * Controlled subjects are defined as explicitly chosen for not having the disease under study
	 * @param chromsWithMut
	 * @param chromsWithoutMut
	 * @param usedChroms 
	 * @return
	 */
	private DataSet assembCaseControlPanel(List<Integer>chromsWithMut, List<Integer>chromsWithoutMut, double snpMaf){
		DataSet uniSet = new DataSet(0, totalCase, totalCtrl);
		Random randomGenerator = new Random();		
		int num_caseAndCarrier = (int)(totalCase*ParamParser.bdr[bdrIndex]);
		//int num_caseNotCarrier = ParamParser.num_ca[ccIndex] - num_caseAndCarrier;
		int num_ctrlAndCarrier = (randomGenerator.nextDouble() < (totalCtrl*ParamParser.bdr[bdrIndex]/ParamParser.grr[grrIndex]))?1:0;
		int num_ctrlNotCarrier = totalCtrl - num_ctrlAndCarrier;		
		//calculate the proportion of homozygotes in cases that are carriers
		double prop_homoInCarrier = snpMaf*snpMaf/(snpMaf*snpMaf + 2*snpMaf*(1-snpMaf));	 	
		int caseCt = 0;
		int controlCt = 0;
		Subject newOne = null;
		File subInfoFile = null;
		if(ParamParser.createSingleSet)
			subInfoFile = new File(location, ParamParser.outPrefix + "singleSetSubInfo.csv");
		else
			subInfoFile = new File(location, ParamParser.outPrefix + "allSubInfo.csv");
		FileHandler handleSubInfoFile = null;
  		if(subInfoFile.exists()){
  			List<Subject> caseCarriers = new ArrayList<Subject>();
  			List<Subject> caseNonCarriers = new ArrayList<Subject>();
  			List<Subject> ctrlCarriers = new ArrayList<Subject>();
  			List<Subject> ctrlNonCarriers = new ArrayList<Subject>();
  			handleSubInfoFile = new FileHandler(subInfoFile);
  			handleSubInfoFile.openForRead();
  			handleSubInfoFile.nextLine();//header
  			String nextLine = null;
  			while((nextLine = handleSubInfoFile.nextLine()) != null){
  				String[] s = nextLine.split(",");
  				int caseOrCtrl = Integer.valueOf(s[3]);
  				int carrierOrNot = Integer.valueOf(s[4]);
  				Subject sub = new Subject(Integer.valueOf(s[0]), Integer.valueOf(s[1]), Integer.valueOf(s[2]), carrierOrNot);
  				if(caseOrCtrl==1 && carrierOrNot==0)
  					caseNonCarriers.add(sub);
  				if(caseOrCtrl==1 && carrierOrNot==1)
  					caseCarriers.add(sub);
  				if(caseOrCtrl==0 && carrierOrNot==0)
  					ctrlNonCarriers.add(sub);
  				if(caseOrCtrl==0 && carrierOrNot==1)
  					ctrlCarriers.add(sub);
  			}
  			handleSubInfoFile.closeReader();
  			
  			int subID = 1;
  			File samplePool = new File(location, ParamParser.outPrefix + "samplePool.csv");
  			if(caseCarriers.size() + caseNonCarriers.size() > totalCase){
  				if(!samplePool.exists()) subInfoFile.renameTo(samplePool);
  				Collections.shuffle(caseCarriers);
  				Collections.shuffle(caseNonCarriers);
  				while(caseCt < num_caseAndCarrier){				
  					newOne = caseCarriers.get(caseCt);
  					newOne.setID(subID++);
  					uniSet.caseList.add(newOne);
  					caseCt++;
  					usedChroms.add(new Integer(newOne.getChromA()));
  					usedChroms.add(new Integer(newOne.getChromB()));
  				}
  				while(caseCt < totalCase){
  					newOne = caseNonCarriers.get(totalCase-caseCt-1);
  					newOne.setID(subID++);
  					uniSet.caseList.add(newOne);
  					caseCt++;
  					usedChroms.add(new Integer(newOne.getChromA()));
  					usedChroms.add(new Integer(newOne.getChromB()));
  				}
  			}
  			else{
  				uniSet.caseList.addAll(caseNonCarriers);
  				uniSet.caseList.addAll(caseCarriers);
  			}
  			
  			if(ctrlCarriers.size() + ctrlNonCarriers.size() > totalCtrl){
  				if(!samplePool.exists()) subInfoFile.renameTo(samplePool);
  				Collections.shuffle(ctrlCarriers);
  				Collections.shuffle(ctrlNonCarriers);
  				while(controlCt < num_ctrlNotCarrier){
  					newOne = ctrlNonCarriers.get(controlCt);
  					newOne.setID(subID++);
  					uniSet.ctrlList.add(newOne);
  					controlCt++;
  					usedChroms.add(new Integer(newOne.getChromA()));
  					usedChroms.add(new Integer(newOne.getChromB()));
  				}
  				while(controlCt < totalCtrl){
  					newOne = ctrlCarriers.get(totalCtrl-controlCt-1);
  					newOne.setID(subID++);
  					uniSet.ctrlList.add(newOne);
  					controlCt++;
  					usedChroms.add(new Integer(newOne.getChromA()));
  					usedChroms.add(new Integer(newOne.getChromB()));
  				}
  			}
  			else{
  				uniSet.ctrlList.addAll(ctrlNonCarriers);
  				uniSet.ctrlList.addAll(ctrlCarriers);
  			}
  			/*if(uniSet.caseList.size() != totalCase || uniSet.ctrlList.size() != totalCtrl){
  				System.out.println("Error: " + uniSet.caseList.size() + " out of " + totalCase + " cases and " + 
  						uniSet.ctrlList.size() + "out of " + totalCtrl + " ctrls are imported.");
  				System.exit(1);
  			}*/
  			Collections.sort(uniSet.caseList);
			Collections.sort(uniSet.ctrlList);
			if(usedChroms.size()==0) importUsedChroms();
			else recordUsedChroms();
			if(samplePool.exists()) samplePool.delete();
  		}
  		else{
			int numOfChromWithMut = chromsWithMut.size();
			int numOfChromWithoutMut = chromsWithoutMut.size();		
			Collections.shuffle(chromsWithMut);
			Collections.shuffle(chromsWithoutMut);
			//int curMutChromIndex = -1;
			//int curWildChromIndex = -1;
			randomGenerator = new Random();
			//select cases from carriers WITH REPLACEMENT (sometimes when the causal variant is very rare there may be only one or two chromosomes having it)	
			int subID = 1;
			while(caseCt < num_caseAndCarrier){
				if(randomGenerator.nextDouble() < prop_homoInCarrier){
					newOne = new Subject(subID++, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
							chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)), 1);
					/*newOne = new Subject("Case" + caseCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
							chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)), 1);*/
					uniSet.caseList.add(newOne);
				}
				else{
					newOne =new Subject(subID++, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
							chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 1);
					/*newOne =new Subject("Case" + caseCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
							chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 1);*/
					uniSet.caseList.add(newOne);	
				}
				caseCt++;
				usedChroms.add(new Integer(newOne.getChromA()));
				usedChroms.add(new Integer(newOne.getChromB()));
			}
			
			//select cases from non-carriers WITH REPLACEMENT
			while(caseCt < totalCase){
				newOne = new Subject(subID++, chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)),
						chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 0);
				/*newOne = new Subject("Case" + caseCt, chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)),
						chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 0);*/
				uniSet.caseList.add(newOne);
				caseCt++;
				usedChroms.add(new Integer(newOne.getChromA()));
				usedChroms.add(new Integer(newOne.getChromB()));
			}
			
			//select controls from non-carriers WITH REPLACEMENT		
			while(controlCt < num_ctrlNotCarrier){
				newOne = new Subject(subID++, chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)),
					chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 0);
				/*newOne = new Subject("Control" + controlCt, chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)),
						chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 0);*/
				uniSet.ctrlList.add(newOne);
				controlCt++;
				usedChroms.add(new Integer(newOne.getChromA()));
				usedChroms.add(new Integer(newOne.getChromB()));
			}
			
			//select controls from carriers WITH REPLACEMENT		
			while(controlCt < totalCtrl){
				if(randomGenerator.nextDouble() < prop_homoInCarrier){
					newOne = new Subject(subID++, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
							chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)), 1);
					/*newOne = new Subject("Control" + controlCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
							chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)), 1);*/
					uniSet.ctrlList.add(newOne);
				}
				else {
					newOne = new Subject(subID++, chromsWithMut.get(Util.randInt(0,numOfChromWithMut-1)),
							chromsWithoutMut.get(Util.randInt(0,numOfChromWithoutMut-1)), 1);
					/*newOne = new Subject("Control" + controlCt, chromsWithMut.get(curMutChromIndex = getNextIndex(chromsWithMut, curMutChromIndex)),
							chromsWithoutMut.get(curWildChromIndex = getNextIndex(chromsWithoutMut, curWildChromIndex)), 1);*/
					uniSet.ctrlList.add(newOne);	
				}
				controlCt++;
				usedChroms.add(new Integer(newOne.getChromA()));
				usedChroms.add(new Integer(newOne.getChromB()));
			}			
			Collections.sort(uniSet.caseList);
			Collections.sort(uniSet.ctrlList);
			recordUsedChroms();
  		}
  		if(!subInfoFile.exists()){
	  		handleSubInfoFile = new FileHandler(subInfoFile);
			handleSubInfoFile.openForWrite(true);
			handleSubInfoFile.writeSubjectFile(uniSet.caseList, uniSet.ctrlList);
			handleSubInfoFile.preventFutureWriting();
			handleSubInfoFile.closeWriter();
  		}
		return uniSet;
	}

	public int getTotalCtrl() {
		return totalCtrl;
	}

	public int getTotalCase() {
		return totalCase;
	}

	public void importUsedChroms() {
		File panelChromsFile = new File(location, "PanelUsedChroms.txt");
		FileHandler panelFileHandler = new FileHandler(panelChromsFile);
		panelFileHandler.openForRead();
		panelFileHandler.nextLine();//header
		String nextLine = null;
		while((nextLine = panelFileHandler.nextLine()) != null)
			usedChroms.add(Integer.valueOf(nextLine));
		panelFileHandler.closeReader();
	}
	
	public void recordUsedChroms(){
	  	File panelChromsFile = new File(location, "PanelUsedChroms.txt");
		if(panelChromsFile.exists()) panelChromsFile.delete();
		FileHandler panelFileHandler = new FileHandler(panelChromsFile);
		panelFileHandler.openForWrite(true);
		panelFileHandler.writeString("IDs of selected chromosomes");
		for(Integer chrom : usedChroms)
			panelFileHandler.writeString(chrom.toString());
		panelFileHandler.closeWriter();
	}

	public Set<Integer> getUsedChroms() {
		return usedChroms;
	}
	
	//For the sake of simplicity, no requirements on the number of homozygotes
/*	public void assembCaseControlPanel(List<Carrier> carriers, List<NonCarrier> nonCarriers, List<char[][]> caseList, List<char[][]> controlList, 
			int num_caseAndCarrier, int num_ctrlNotCarrier){		
		int caseCt = 0;
		int controlCt = 0;	
		Collections.shuffle(carriers);
		Collections.shuffle(nonCarriers);
		int curCarrierIndex = 0;
		int curNonCarrierIndex = 0;
		while(caseCt < num_caseAndCarrier){
			if(curCarrierIndex == carriers.size()){
				Collections.shuffle(carriers);
				curCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];//must be a different object
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(0), 1, diploid[0], 0, ParamParser.num_mark);
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(1), 1, diploid[1], 0, ParamParser.num_mark);
			caseList.add(diploid);
			curCarrierIndex++;
			caseCt++;
		}
		while(caseCt < ParamParser.num_ca[caseIndex]){
			if(curNonCarrierIndex == nonCarriers.size()){
				Collections.shuffle(nonCarriers);
				curNonCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];
			diploid[0] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(0).clone();
			diploid[1] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(1).clone();
			caseList.add(diploid);
			curNonCarrierIndex++;
			caseCt++;
		}				
		while(controlCt < num_ctrlNotCarrier){
			if(curNonCarrierIndex == nonCarriers.size()){
				Collections.shuffle(nonCarriers);
				curNonCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];
			diploid[0] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(0).clone();
			diploid[1] = nonCarriers.get(curNonCarrierIndex).getSingleChrom(1).clone();
			controlList.add(diploid);
			curNonCarrierIndex++;
			controlCt++;
		}				
		while(controlCt < ParamParser.num_co[ctrlIndex]){
			if(curCarrierIndex == carriers.size()){
				Collections.shuffle(carriers);
				curCarrierIndex = 0;
			}
			char[][] diploid = new char[2][ParamParser.num_mark];
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(0), 1, diploid[0], 0, ParamParser.num_mark);
			System.arraycopy(carriers.get(curCarrierIndex).getSingleChrom(1), 1, diploid[1], 0, ParamParser.num_mark);
			controlList.add(diploid);
			curCarrierIndex++;
			controlCt++;
		}		
	}*/
}
