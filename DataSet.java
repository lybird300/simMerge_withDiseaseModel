package seqSIM;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;

public class DataSet {
	public List<Subject> caseList;
	public List<Subject> ctrlList;
	public List<Integer> tagList;
	public int dataSetID;
	public Set<Integer> snpWithAlleleSwitched = null;
	public Set<Integer> snpWithGenoShuffled = null;
	//public List<Integer> missingGeno = null;
	//public List<Integer> wrongGeno = null;
	public Map<Integer, int[]> missingGenoList = null;//<absolute ID of missing genotype, [subjectID, tag 20M position]>
	public Map<Integer, int[]> wrongGenoList = null;//<absolute ID of missing genotype, [subjectID, tag 20M position, correct geno, wrong geno]>
	
	public DataSet(int dataSetID, int num_case, int num_ctrl){
		this.dataSetID = dataSetID;
		this.caseList = new ArrayList<Subject>(num_case);
		this.ctrlList = new ArrayList<Subject>(num_ctrl);
	}
	
	public int getTagCt(){
		return tagList.size();
	}
	
	public int getDataSetSize(){
		return caseList.size() + ctrlList.size();
	}
	
	public List<Integer> getTagIDs(){
		List<Integer> tagIDList = new ArrayList<Integer>();
		for(int i = 0; i < tagList.size(); i++){
			tagIDList.add(i);
		}
		return tagIDList;
	}

	public List<Integer> getTagList() {
		return tagList;
	}

	public void assignErr(int errType, int sizeOfPool, Set<Integer> itemToExclude, int numOfItemsToSelect) {
		Set<Integer> selectedItems = new HashSet<Integer>();
		Random randGenerator = new Random();
		while(selectedItems.size() < numOfItemsToSelect){
			int i = randGenerator.nextInt(sizeOfPool);
			if(itemToExclude!=null && itemToExclude.contains(i)) continue;
			else selectedItems.add(i);
		}
		List<Integer> sortedItems = new ArrayList<Integer>(selectedItems);
		selectedItems.clear();
		selectedItems = null;
		Collections.sort(sortedItems);
		switch(errType){
			case 0:	snpWithAlleleSwitched = new HashSet<Integer>(sortedItems);
					break;
			case 1:	snpWithGenoShuffled = new HashSet<Integer>(sortedItems);
					break;
			case 2: missingGenoList = new HashMap<Integer, int[]>(numOfItemsToSelect);
					for(int item : sortedItems)
						missingGenoList.put(item, null);
					break;
			case 3: wrongGenoList = new HashMap<Integer, int[]>(numOfItemsToSelect);
					for(int item : sortedItems)
						wrongGenoList.put(item, null);
					break;
		}
	}

	public Map<Integer, int[]> getMissingGeno() {
		return missingGenoList;
	}

	public Map<Integer, int[]> getWrongGeno() {
		return wrongGenoList;
	}
	
/*	public void outputLgenFile(File lgenFile, byte[][] haplotypes, Map<Integer, Integer> tag20mPosTo2MIndex, List<Map<Character, Character>> recodeMap) {
		FileHandler writeLgen = new FileHandler(lgenFile);
		writeLgen.openForWrite(false);
		int subUniID;
		int subID = 0;
		int ATsum = 0;
		int CGsum = 0;
		for(;subID<caseList.size();subID++){
			subUniID = caseList.get(subID).getID();
			for(int snpID = 0; snpID < tagList.size(); snpID++){
				int tag20MPos = tagList.get(snpID);
				Character[] alleles = recode(tag20MPos, haplotypes[2*subID][snpID], haplotypes[2*subID+1][snpID], tag20mPosTo2MIndex, recodeMap);
				if((alleles[0]=='A' && alleles[1]=='T')|| (alleles[0]=='T' && alleles[1]=='A')) ATsum++;
				if((alleles[0]=='C' && alleles[1]=='G')|| (alleles[0]=='G' && alleles[1]=='C')) CGsum++;
				writeLgen.writeString(subUniID + "\t" + subUniID + "\t" + "rs" + tag20MPos + "\t" + alleles[0] + "\t" + alleles[1]);
			}
		}
		for(int ct = 0;ct<ctrlList.size();ct++){
			subUniID = ctrlList.get(ct).getID();
			for(int snpID = 0; snpID < tagList.size(); snpID++){
				int tag20MPos = tagList.get(snpID);
				Character[] alleles = recode(tag20MPos, haplotypes[2*(subID+ct)][snpID], haplotypes[2*(subID+ct)+1][snpID], tag20mPosTo2MIndex, recodeMap);
				if((alleles[0]=='A' && alleles[1]=='T')|| (alleles[0]=='T' && alleles[1]=='A')) ATsum++;
				if((alleles[0]=='C' && alleles[1]=='G')|| (alleles[0]=='G' && alleles[1]=='C')) CGsum++;
				if(ct==ctrlList.size()-1 && snpID == tagList.size()-1)
					writeLgen.writeLastLine(subUniID + "\t" + subUniID + "\t" + "rs" + tag20MPos + "\t" + alleles[0] + "\t" + alleles[1]);
				else 
					writeLgen.writeString(subUniID + "\t" + subUniID + "\t" + "rs" + tag20MPos + "\t" + alleles[0] + "\t" + alleles[1]);
			}
		}
		System.out.println("There are " + ATsum + " AT genotypes and " + CGsum + " CG genotypes.\n"
				+ "The proportion of AT/CG genotypes are " + (ATsum+CGsum)/(double)tagList.size());
		writeLgen.closeWriter();
	}
	
	private Character[] recode(int tag20mPos, byte oldmAllele, byte oldpAllele, Map<Integer, Integer> tag20mPosTo2MIndex,
			List<Map<Character, Character>> recodeMap) {
		Character mAllele = Byte.toString(oldmAllele).charAt(0);
		Character pAllele = Byte.toString(oldpAllele).charAt(0);
		Character[] alleles = new Character[2];
		if(recodeMap==null){
			alleles[0] = mAllele;
			alleles[1] = pAllele;	
		}
		else{
			alleles[0] = recodeMap.get(tag20mPosTo2MIndex.get(tag20mPos)).get(mAllele);
			alleles[1] = recodeMap.get(tag20mPosTo2MIndex.get(tag20mPos)).get(pAllele);		
		}
		return alleles;
	}*/
	
	private Character[][] recode(byte[][] haplotypes, Map<Integer, Integer> tag20mPosTo2MIndex, List<Map<Character, Character>> recodeMap){
		Character[][] recodedHaps = new Character[haplotypes.length][haplotypes[0].length];
		for(int snpID = 0; snpID < tagList.size(); snpID++){
			Map<Character, Character> snpMap = recodeMap.get(tag20mPosTo2MIndex.get(tagList.get(snpID)));
			for(int hapCt = 0; hapCt < haplotypes.length; hapCt++)
				recodedHaps[hapCt][snpID] = snpMap.get(Byte.toString(haplotypes[hapCt][snpID]).charAt(0));
		}
		return recodedHaps;
	}
	
	private void outputLgenFile(File lgenFile, byte[][] haplotypes) {
		FileHandler writeLgen = new FileHandler(lgenFile);
		writeLgen.openForWrite(false);
		int subUniID;
		int subID = 0;
		for(;subID<caseList.size();subID++){
			subUniID = caseList.get(subID).getID();
			for(int snpID = 0; snpID < tagList.size(); snpID++)
				writeLgen.writeString(subUniID + "\t" + subUniID + "\t" + "rs" + tagList.get(snpID) + "\t" + String.valueOf(haplotypes[2*subID][snpID]) + 
						"\t" + String.valueOf(haplotypes[2*subID+1][snpID]));
		}
		for(int ct = 0;ct<ctrlList.size();ct++){
			subUniID = ctrlList.get(ct).getID();
			for(int snpID = 0; snpID < tagList.size(); snpID++){
				if(ct==ctrlList.size()-1 && snpID == tagList.size()-1)
					writeLgen.writeLastLine(subUniID + "\t" + subUniID + "\t" + "rs" + tagList.get(snpID) + "\t" + String.valueOf(haplotypes[2*(subID+ct)][snpID]) + 
							"\t" + String.valueOf(haplotypes[2*(subID+ct)+1][snpID]));
				else 
					writeLgen.writeString(subUniID + "\t" + subUniID + "\t" + "rs" + tagList.get(snpID) + "\t" + String.valueOf(haplotypes[2*(subID+ct)][snpID]) + 
						"\t" + String.valueOf(haplotypes[2*(subID+ct)+1][snpID]));
			}
		}
		writeLgen.closeWriter();
	}
	
	private void outputLgenFile(File lgenFile, Character[][] haplotypes) {
		FileHandler writeLgen = new FileHandler(lgenFile);
		writeLgen.openForWrite(false);
		int subUniID;
		int subID = 0;
		for(;subID<caseList.size();subID++){
			subUniID = caseList.get(subID).getID();
			for(int snpID = 0; snpID < tagList.size(); snpID++)
				writeLgen.writeString(subUniID + "\t" + subUniID + "\t" + "rs" + tagList.get(snpID) + "\t" + haplotypes[2*subID][snpID] + 
						"\t" + haplotypes[2*subID+1][snpID]);
		}
		for(int ct = 0;ct<ctrlList.size();ct++){
			subUniID = ctrlList.get(ct).getID();
			for(int snpID = 0; snpID < tagList.size(); snpID++){
				if(ct==ctrlList.size()-1 && snpID == tagList.size()-1)
					writeLgen.writeLastLine(subUniID + "\t" + subUniID + "\t" + "rs" + tagList.get(snpID) + "\t" + haplotypes[2*(subID+ct)][snpID] + 
							"\t" + haplotypes[2*(subID+ct)+1][snpID]);
				else 
					writeLgen.writeString(subUniID + "\t" + subUniID + "\t" + "rs" + tagList.get(snpID) + "\t" + haplotypes[2*(subID+ct)][snpID] + 
						"\t" + haplotypes[2*(subID+ct)+1][snpID]);
			}
		}
		writeLgen.closeWriter();
	}

	public byte[][] getUniSetCorrect12Haplotypes(String workDir, Map<Integer, Integer> tagPosToIndex, Integer[] usedChromsArray, 
			Map<Integer, Integer> tag20MPosTo2MIndex, List<Map<Character, Character>> recodeMap, int[][] subIDtoHapID){
		Integer[] tag2MIndex = tag20MPosTo2MIndex.values().toArray(new Integer[tag20MPosTo2MIndex.size()]);//This is based on the unified set of tags
		File trueLgen = new File(workDir, ParamParser.outPrefix + "Set" + dataSetID + "_noErr12.lgen");
		byte[][] haplotypes = null;
		if(trueLgen.exists())
			haplotypes = convertLgenFileToHapMatrix(trueLgen);
		else{
			Map<Integer, byte[]> chromIDToSeq = Util.gatherMutMarkersfromDistFiles(usedChromsArray, tag2MIndex);
			haplotypes = new byte[2*this.getDataSetSize()][this.getTagCt()];
			byte[] mHap = null;
			byte[] pHap = null;
			for(Subject csSub : caseList){
				int cellID = csSub.getID() - 1;
				mHap = chromIDToSeq.get(subIDtoHapID[cellID][0]);
				pHap = chromIDToSeq.get(subIDtoHapID[cellID][1]);
				for(int snpID = 0; snpID < tagList.size(); snpID++){
					haplotypes[cellID*2][snpID] =  mHap[tagPosToIndex.get(tagList.get(snpID))];
					haplotypes[cellID*2+1][snpID] =  pHap[tagPosToIndex.get(tagList.get(snpID))];
				}
			}
			for(Subject ctSub : ctrlList){
				int cellID = ctSub.getID() - 1;
				mHap = chromIDToSeq.get(subIDtoHapID[cellID][0]);
				pHap = chromIDToSeq.get(subIDtoHapID[cellID][1]);
				for(int snpID = 0; snpID < tagList.size(); snpID++){
					haplotypes[cellID*2][snpID] =  mHap[tagPosToIndex.get(tagList.get(snpID))];
					haplotypes[cellID*2+1][snpID] =  pHap[tagPosToIndex.get(tagList.get(snpID))];
				}
			}		
			outputLgenFile(trueLgen, haplotypes);
		}
		File recodedLgen = new File(workDir, ParamParser.outPrefix + "Set" + dataSetID + "_noErr.lgen");
		if(!recodedLgen.exists())
			outputLgenFile(recodedLgen, recode(haplotypes, tag20MPosTo2MIndex, recodeMap));
		return haplotypes;		
	}
	
	private byte[][] convertLgenFileToHapMatrix(File trueLgen) {
		byte[][] haplotypes = new byte[2*this.getDataSetSize()][this.getTagCt()];
		FileHandler handleLgenFile = new FileHandler(trueLgen);
		handleLgenFile.openForRead();
		String nextLine = null;
		int markerID = 0;
		int subjectSubSetID = 0;
		while((nextLine = handleLgenFile.nextLine()) != null){
			String[] s = nextLine.split("\t");		
			haplotypes[subjectSubSetID*2][markerID] = Byte.valueOf(s[3]);
			haplotypes[subjectSubSetID*2+1][markerID] = Byte.valueOf(s[4]);
			if((++markerID)>=haplotypes[0].length){
				subjectSubSetID++;
				markerID = 0;
			}
		}
		handleLgenFile.closeReader();
		return haplotypes;
	}

	public void generateSubSetErrHaplotypes(String workDir, byte[][] subsetCorrectHaplotypes, Map<Integer, Integer> tag20mPosTo2MIndex, 
			List<Map<Character, Character>> recodeMap, ErrController ec){
		File lgenWithErr = new File(workDir, ParamParser.outPrefix + "Set" + dataSetID + "_flipErr12.lgen");
		byte[][] haplotypesWithErr = null;
		if(lgenWithErr.exists())
			haplotypesWithErr = convertLgenFileToHapMatrix(lgenWithErr);
		else{
			haplotypesWithErr = new byte[subsetCorrectHaplotypes.length][subsetCorrectHaplotypes[0].length];
			Map<Integer, Byte> snpWithMinorHomoChangeToHetero = new HashMap<Integer, Byte>();
			if(ParamParser.randMarkerRate > 0 && ParamParser.focusOnMinorHomo && snpWithGenoShuffled==null){
				Map<Integer, Byte> markersWithMinorHomoSubjects = new HashMap<Integer, Byte>();
				for(int tagID = 0; tagID < this.getTagCt(); tagID++){
					if(ParamParser.imposeErrOnATCGMarkerOnly && !isATCGMarker(recodeMap.get(tag20mPosTo2MIndex.get(tagList.get(tagID))))) continue;
					if(!ParamParser.imposeErrOnATCGMarkerOnly && isATCGMarker(recodeMap.get(tag20mPosTo2MIndex.get(tagList.get(tagID))))) continue;
					Set<Integer> homoGroup1 = new HashSet<Integer>();
					Set<Integer> homoGroup2 = new HashSet<Integer>();
					for(int s = 0; s < this.getDataSetSize(); s++){
						if((subsetCorrectHaplotypes[2*s][tagID]== subsetCorrectHaplotypes[2*s+1][tagID])){
							if(subsetCorrectHaplotypes[2*s][tagID]==1) homoGroup1.add(s);
							else homoGroup2.add(s);
						}
					}
					if(homoGroup1.size() < 8 || homoGroup2.size() < 8) continue;
					Byte minorAllele = (homoGroup1.size()< homoGroup2.size())? (byte)1:(byte)2;
					markersWithMinorHomoSubjects.put(tagID, minorAllele);
				}
				List<Integer> candidateSNPs = new ArrayList<Integer>(markersWithMinorHomoSubjects.keySet());
				int numOfGenoShuffledSNPs = (int)(ParamParser.randMarkerRate*this.getTagCt());
				snpWithGenoShuffled = new HashSet<Integer>();
				Collections.shuffle(candidateSNPs);
				for(int i = 0; i < numOfGenoShuffledSNPs; i++){
					Integer oneSNP = candidateSNPs.get(i);
					snpWithGenoShuffled.add(oneSNP);
					snpWithMinorHomoChangeToHetero.put(oneSNP, markersWithMinorHomoSubjects.get(oneSNP));
				}
				File shuffledFile = new File(workDir, ParamParser.outPrefix + "ShuffledSNPs_Set" + dataSetID + ".txt");
				ec.allShuffledSNPs.addAll(recordGenoShuffledSNPs(shuffledFile));
			}
			for(int j = 0; j < this.getTagCt(); j++){
				int genoCt = 0;
				boolean doGenoShuffle = false;
				Map<Integer, Integer> oldSubToNewSub = null;
				Set<Integer> selectedSubjects = null;
				if(ParamParser.randMarkerRate > 0 && snpWithGenoShuffled.contains(j)){
					//System.out.println("Set" + dataSetID + " shuffled marker pos: " + tagList.get(j));
					if(ParamParser.completeRandShuffle){
						oldSubToNewSub = new HashMap<Integer, Integer>();
						List<Integer> genoPerSNPIDMask = new ArrayList<Integer>(getDataSetSize());
						for(int g = 0; g < getDataSetSize(); g++)
							genoPerSNPIDMask.add(g);
						Collections.shuffle(genoPerSNPIDMask);
						for(int g = 0; g < getDataSetSize(); g++){
							if(!oldSubToNewSub.containsKey(g)) oldSubToNewSub.put(g, genoPerSNPIDMask.get(g));
						}
					}
					else if(ParamParser.startFromHomozygote){
						if(ParamParser.focusOnMinorHomo){
							List<Integer> minorHomoID = new ArrayList<Integer>();
							for(int s = 0; s < this.getDataSetSize(); s++){
								if((subsetCorrectHaplotypes[2*s][j]== subsetCorrectHaplotypes[2*s+1][j]) && 
										subsetCorrectHaplotypes[2*s][j]==snpWithMinorHomoChangeToHetero.get(j)){
									minorHomoID.add(s);
								}
							}
							int numOfSubjectsNeedToBeChanged = (int)(ParamParser.subSwapRate*minorHomoID.size());
							selectedSubjects = new HashSet<Integer>();
							Collections.shuffle(minorHomoID);
							for(int i = 0; i < numOfSubjectsNeedToBeChanged; i++){
								selectedSubjects.add(minorHomoID.get(i));
							}
						}else{
							oldSubToNewSub = new HashMap<Integer, Integer>();
							List<Integer> homozygoteID = new ArrayList<Integer>();
							List<Integer> heterozygoteID = new ArrayList<Integer>();
							for(int s = 0; s < this.getDataSetSize(); s++){
								if(subsetCorrectHaplotypes[2*s][j]== subsetCorrectHaplotypes[2*s+1][j])
									homozygoteID.add(s);
								else
									heterozygoteID.add(s);
							}
							int numOfSubjectsNeedToBeSwapped = (int)(ParamParser.subSwapRate*getDataSetSize());
							int numOfSubjectsSwapped = 0; 
							while(numOfSubjectsSwapped < numOfSubjectsNeedToBeSwapped){
								Collections.shuffle(homozygoteID);
								Collections.shuffle(heterozygoteID);
								int subID1 = homozygoteID.get(0);
								int subID2 = heterozygoteID.get(0);
								if(!oldSubToNewSub.containsKey(subID1) && !oldSubToNewSub.containsKey(subID2)){
									oldSubToNewSub.put(subID1, subID2);
									oldSubToNewSub.put(subID2, subID1);
									numOfSubjectsSwapped += 2;
								}
								else if(ParamParser.allowDoubleSwapErr){
									if(!oldSubToNewSub.containsKey(subID1) || !oldSubToNewSub.containsKey(subID2))
										numOfSubjectsSwapped++;
									int newIDForSubID1 = (oldSubToNewSub.containsKey(subID2))? (oldSubToNewSub.get(subID2)):subID2;
									int newIDForSubID2 = (oldSubToNewSub.containsKey(subID1))? (oldSubToNewSub.get(subID1)):subID1;
									oldSubToNewSub.put(subID1, newIDForSubID1);
									oldSubToNewSub.put(subID2, newIDForSubID2);
								}
							}
						}
					}else{
						int numOfSubjectsNeedToBeChanged = (int)(ParamParser.subSwapRate*getDataSetSize());
						int subjectCt = this.getDataSetSize();
						selectedSubjects = new HashSet<Integer>();
						while(selectedSubjects.size() < numOfSubjectsNeedToBeChanged){
							int i = (new Random()).nextInt(subjectCt);
							selectedSubjects.add(i);
						}
					}
					doGenoShuffle = true;
				}
				for(int i = 0; i < this.getDataSetSize(); i++, genoCt++){
					if(ParamParser.missingRate > 0 && missingGenoList.containsKey(genoCt)){
						haplotypesWithErr[2*i][j] = (byte)0;
						haplotypesWithErr[2*i+1][j] = (byte)0;
			    		int[] missingGenoRecord = missingGenoList.get(genoCt);
			    		if(missingGenoRecord == null){
			    			missingGenoRecord = new int[2];
				    		int subUniID = (i >= caseList.size())?ctrlList.get(i-caseList.size()).getID():caseList.get(i).getID();
				    		missingGenoRecord[0] = subUniID;
				    		missingGenoRecord[1] = tagList.get(j);
				    		missingGenoList.put(genoCt, missingGenoRecord);
			    		}
						continue;
					}
					if(ParamParser.switchRate > 0){
						if(!snpWithAlleleSwitched.contains(j)){
							haplotypesWithErr[2*i][j] = subsetCorrectHaplotypes[2*i][j];
							haplotypesWithErr[2*i+1][j] = subsetCorrectHaplotypes[2*i+1][j];
						}
						else if((ParamParser.imposeErrBasedOAlleleStates && ParamParser.imposeErrOnATCGMarkerOnly) || 
							(ParamParser.imposeErrBasedOAlleleStates && !ParamParser.imposeErrOnATCGMarkerOnly)){
							haplotypesWithErr[2*i][j] = (byte)(3-subsetCorrectHaplotypes[2*i][j]);
							haplotypesWithErr[2*i+1][j] = (byte)(3-subsetCorrectHaplotypes[2*i+1][j]);
						}
					}
					if(ParamParser.randMarkerRate > 0){
						if(doGenoShuffle){
							if(oldSubToNewSub!=null && selectedSubjects==null){
								int newi = (oldSubToNewSub.containsKey(i))?(oldSubToNewSub.get(i)):i;
								haplotypesWithErr[2*i][j] = subsetCorrectHaplotypes[2*newi][j];
								haplotypesWithErr[2*i+1][j] = subsetCorrectHaplotypes[2*newi+1][j];
							}
							if(oldSubToNewSub==null && selectedSubjects!=null){
								if(selectedSubjects.contains(i)){
									if (!ParamParser.focusOnMinorHomo) {
										if (subsetCorrectHaplotypes[2 * i][j] == subsetCorrectHaplotypes[2 * i + 1][j]) {
											haplotypesWithErr[2 * i][j] = subsetCorrectHaplotypes[2 * i][j];
											haplotypesWithErr[2 * i + 1][j] = (byte) (3 - subsetCorrectHaplotypes[2 * i + 1][j]);
										} else {
											Random r = new Random();
											Byte allele = (r.nextDouble() < 0.5) ? (byte) 1 : (byte) 2;
											haplotypesWithErr[2 * i][j] = allele.byteValue();
											haplotypesWithErr[2 * i + 1][j] = allele.byteValue();
										} 
									}
									else{
										haplotypesWithErr[2 * i][j] = subsetCorrectHaplotypes[2 * i][j];
										haplotypesWithErr[2 * i + 1][j] = (byte) (3 - subsetCorrectHaplotypes[2 * i + 1][j]);
									}
								}else{
									haplotypesWithErr[2*i][j] = subsetCorrectHaplotypes[2*i][j];
									haplotypesWithErr[2*i+1][j] = subsetCorrectHaplotypes[2*i+1][j];
								}
							}
						}else{
							haplotypesWithErr[2*i][j] = subsetCorrectHaplotypes[2*i][j];
							haplotypesWithErr[2*i+1][j] = subsetCorrectHaplotypes[2*i+1][j];
						}
					}
					if(ParamParser.errorRate > 0 && wrongGenoList.containsKey(genoCt)){					
						int[] wrongGenoRecord = wrongGenoList.get(genoCt);
						if(wrongGenoRecord == null){
							wrongGenoRecord = new int[6];
							int subUniID = (i >= caseList.size())?ctrlList.get(i-caseList.size()).getID():caseList.get(i).getID();
							wrongGenoRecord[0] = subUniID;
							wrongGenoRecord[1] = tagList.get(j);
							wrongGenoRecord[2] = haplotypesWithErr[2*i][j];
							wrongGenoRecord[3] = haplotypesWithErr[2*i+1][j];
							byte[] changedAlleles = getWrongGenoBasedOnAlleleFreq(subsetCorrectHaplotypes[2*i][j], subsetCorrectHaplotypes[2*i+1][j], j, subsetCorrectHaplotypes);
							haplotypesWithErr[2*i][j] = changedAlleles[0];
							haplotypesWithErr[2*i+1][j] = changedAlleles[1];
							wrongGenoRecord[4] = changedAlleles[0];
							wrongGenoRecord[5] = changedAlleles[1];
							wrongGenoList.put(genoCt, wrongGenoRecord);
						}
						else{
							haplotypesWithErr[2*i][j] = (byte)wrongGenoRecord[4];
							haplotypesWithErr[2*i+1][j] = (byte)wrongGenoRecord[5];
						}
					}
				}
			}
			this.outputLgenFile(lgenWithErr, haplotypesWithErr);
		}
    	File recodedErrLgen = new File(workDir, ParamParser.outPrefix + "Set" + dataSetID + "_flipErr.lgen");
    	if(!recodedErrLgen.exists())
    		outputLgenFile(recodedErrLgen, recode(haplotypesWithErr, tag20mPosTo2MIndex, recodeMap));	
	}

	private byte[] getWrongGenoBasedOnAlleleFreq(byte originMAllele, byte originPAllele, int snpID, byte[][] correctSubSetHaplotypes){
	    byte[] changedAlleles = {originMAllele, originPAllele};	    
	    int originMAlleleCt = 0;
	    for (int i = 0; i < correctSubSetHaplotypes.length; i++){
	      if(correctSubSetHaplotypes[i][snpID] == originMAllele)
	        originMAlleleCt++;
	    }
		int majorAllele;
		double majorAlleleFreq;
		if(originMAlleleCt >= this.getDataSetSize()){
			majorAllele = originMAllele;
			majorAlleleFreq = originMAlleleCt/(correctSubSetHaplotypes.length*2.0);
		}
		else{
			majorAllele = 3-originMAllele;//1->2, 2->1
			majorAlleleFreq = 1-originMAlleleCt/(correctSubSetHaplotypes.length*2.0);
		}
		Random randGenerator = new Random();
		if(originMAllele==originPAllele && originMAllele==(byte)majorAllele){//current genotype is major allele homozygote
			if(randGenerator.nextDouble() < Math.pow((1-majorAlleleFreq),2)/(Math.pow((1-majorAlleleFreq),2) + 2*majorAlleleFreq*(1-majorAlleleFreq)))
				changedAlleles[0]=changedAlleles[1]=(byte)(3-majorAllele);
			else{
				if(randGenerator.nextDouble() < 0.5) changedAlleles[0] = (byte)(3-majorAllele);
				else changedAlleles[1] = (byte)(3-majorAllele);
			}
		}
		else if(originMAllele==originPAllele && originMAllele!=(byte)majorAllele){//current genotype is minor allele homozygote
			if(randGenerator.nextDouble() > 2*majorAlleleFreq*(1-majorAlleleFreq)/(Math.pow(majorAlleleFreq,2) + 2*majorAlleleFreq*(1-majorAlleleFreq)))
				changedAlleles[0]=changedAlleles[1]=(byte)majorAllele;
			else{
				if(randGenerator.nextDouble() < 0.5) changedAlleles[0] = (byte)majorAllele;
				else changedAlleles[1] = (byte)majorAllele;
			}
		}
		else{//current genotype is heterozygote
			if(randGenerator.nextDouble() > Math.pow(majorAlleleFreq,2)/(Math.pow(majorAlleleFreq,2) + Math.pow((1-majorAlleleFreq),2)))
				changedAlleles[0]=changedAlleles[1]=(byte)(3-majorAllele);
			else
				changedAlleles[0]=changedAlleles[1]=(byte)majorAllele;
		}
		return changedAlleles;
	}

	public void generateUniSetErrHaplotypes(String workDir, Map<Integer, Integer> tagPosToUniTagListIndex, byte[][] uniSetCorrectHaplotypes, 
			ErrController ec, Map<Integer, Integer> tag20MPosTo2MIndex, List<Map<Character, Character>> recodeMap) {
		Integer SetID = 0;
		if(ParamParser.imposeFlipErrOnUniSet && !ParamParser.createSingleSet) SetID = 0;
		if(ParamParser.createSingleSet && !ParamParser.imposeFlipErrOnUniSet) SetID = 1;
	    File lgenWithErr  = new File(workDir, ParamParser.outPrefix + "Set" + SetID.toString() + "_genoShuffleErr12.lgen");
	    byte[][] haplotypesWithErr = null;
		if(lgenWithErr.exists())
			haplotypesWithErr = convertLgenFileToHapMatrix(lgenWithErr);
		else{
			if(ec.allShuffledSNPs==null){
				//TODO: read from record files
				
			}
			haplotypesWithErr = new byte[uniSetCorrectHaplotypes.length][uniSetCorrectHaplotypes[0].length];
			for(int j = 0; j < this.getTagCt(); j++){
				boolean doGenoShuffle = false;
				Map<Integer, Integer> oldSubToNewSub = new HashMap<Integer, Integer>();
				if(ParamParser.randMarkerRate > 0 && ec.allShuffledSNPs.contains(j)){
					if(ParamParser.completeRandShuffle){
						List<Integer> genoPerSNPIDMask = new ArrayList<Integer>(getDataSetSize());
						for(int g = 0; g < getDataSetSize(); g++)
							genoPerSNPIDMask.add(g);
						Collections.shuffle(genoPerSNPIDMask);
						for(int g = 0; g < getDataSetSize(); g++){
							if(!oldSubToNewSub.containsKey(g)) oldSubToNewSub.put(g, genoPerSNPIDMask.get(g));
						}
					}
					else{
						List<Integer> homozygoteID = new ArrayList<Integer>();
						List<Integer> heterozygoteID = new ArrayList<Integer>();
						for(int s = 0; s < this.getDataSetSize(); s++){
							if(uniSetCorrectHaplotypes[2*s][j]== uniSetCorrectHaplotypes[2*s+1][j])
								homozygoteID.add(s);
							else
								heterozygoteID.add(s);
						}
						int numOfSubjectsNeedToBeSwapped = (int)(ParamParser.subSwapRate*getDataSetSize());
						int numOfSubjectsSwapped = 0; 
						while(numOfSubjectsSwapped < numOfSubjectsNeedToBeSwapped){
							Collections.shuffle(homozygoteID);
							Collections.shuffle(heterozygoteID);
							int subID1 = homozygoteID.get(0);
							int subID2 = heterozygoteID.get(0);
							if(!oldSubToNewSub.containsKey(subID1) && !oldSubToNewSub.containsKey(subID2)){
								oldSubToNewSub.put(subID1, subID2);
								oldSubToNewSub.put(subID2, subID1);
								numOfSubjectsSwapped += 2;
							}
							else if(ParamParser.allowDoubleSwapErr){
								if(!oldSubToNewSub.containsKey(subID1) || !oldSubToNewSub.containsKey(subID2))
									numOfSubjectsSwapped++;
								int newIDForSubID1 = (oldSubToNewSub.containsKey(subID2))? (oldSubToNewSub.get(subID2)):subID2;
								int newIDForSubID2 = (oldSubToNewSub.containsKey(subID1))? (oldSubToNewSub.get(subID1)):subID1;
								oldSubToNewSub.put(subID1, newIDForSubID1);
								oldSubToNewSub.put(subID2, newIDForSubID2);
							}
						}
					}
					System.out.println("Actually swapped genotypes: " + oldSubToNewSub.size());
					doGenoShuffle = true;
				}
				for(int i = 0; i < this.getDataSetSize(); i++){
					if(!ParamParser.imposeFlipErrOnUniSet && !ParamParser.createSingleSet){
						haplotypesWithErr[2*i][j] = uniSetCorrectHaplotypes[2*i][j];
						haplotypesWithErr[2*i+1][j] = uniSetCorrectHaplotypes[2*i+1][j];
					}
					else{
						if(snpWithAlleleSwitched.contains(j)){
							haplotypesWithErr[2*i][j] = (byte)(3-uniSetCorrectHaplotypes[2*i][j]);
							haplotypesWithErr[2*i+1][j] = (byte)(3-uniSetCorrectHaplotypes[2*i+1][j]);
						}
						else{
							haplotypesWithErr[2*i][j] = uniSetCorrectHaplotypes[2*i][j];
							haplotypesWithErr[2*i+1][j] = uniSetCorrectHaplotypes[2*i+1][j];
						}
					}
					if(doGenoShuffle){
						int newi = (oldSubToNewSub.containsKey(i))?(oldSubToNewSub.get(i)):i;
						haplotypesWithErr[2*i][j] = uniSetCorrectHaplotypes[2*newi][j];
						haplotypesWithErr[2*i+1][j] = uniSetCorrectHaplotypes[2*newi+1][j];		
					}	
				}
			}
			if(ParamParser.missingRate > 0){
				for(Entry<Integer, int[]>  e : ec.allMissingGeno.entrySet()){
					int[] record = e.getValue();
					int cellID = record[0] - 1;
					int tagID = tagPosToUniTagListIndex.get(record[1]);
					haplotypesWithErr[2*cellID][tagID] = 0;
					haplotypesWithErr[2*cellID+1][tagID] = 0;
				}
			}
			if(ParamParser.errorRate > 0){
				for(Entry<Integer, int[]>  e : ec.allWrongGeno.entrySet()){
					int[] record = e.getValue();
					int cellID = record[0] - 1;
					int tagID = tagPosToUniTagListIndex.get(record[1]);
					haplotypesWithErr[2*cellID][tagID] = (byte)record[2];
					haplotypesWithErr[2*cellID+1][tagID] = (byte)record[3];
				}
			}
			this.outputLgenFile(lgenWithErr, haplotypesWithErr);
		}
		File recodedLgenFile = new File(workDir, ParamParser.outPrefix + "Set" + SetID.toString() + "_genoShuffleErr.lgen");//If flip error is the only type of errors, this file is the same as Set0_noErr.lgen
	    if(!recodedLgenFile.exists()) 
	    	this.outputLgenFile(recodedLgenFile, recode(haplotypesWithErr, tag20MPosTo2MIndex, recodeMap));
	}

	public byte[][] getSubSetCorrect12Haplotypes(String workDir, Map<Integer, Integer> tagPosToUniTagListIndex, byte[][] uniSetCorrectGenotype,
			Map<Integer, Integer> tag20MPosTo2MIndex, List<Map<Character, Character>> recodeMap) {
		File trueLgen = new File(workDir, ParamParser.outPrefix + "Set" + dataSetID + "_noErr12.lgen");
		byte[][] haplotypes = null;
		if(trueLgen.exists())
			haplotypes = convertLgenFileToHapMatrix(trueLgen);
		else{
			haplotypes = new byte[2*this.getDataSetSize()][this.getTagCt()];
			for(int j = 0; j < haplotypes[0].length; j++){
				int snpUniID = tagPosToUniTagListIndex.get(tagList.get(j));
				int k = 0;
				for(Subject csSub : caseList){
					int cellID = csSub.getID()-1;
					haplotypes[k*2][j] =  uniSetCorrectGenotype[cellID*2][snpUniID];
					haplotypes[k*2+1][j] =  uniSetCorrectGenotype[cellID*2+1][snpUniID];
					k++;
				}
				for(Subject ctSub : ctrlList){
					int cellID = ctSub.getID() - 1;
					haplotypes[k*2][j] =  uniSetCorrectGenotype[cellID*2][snpUniID];
					haplotypes[k*2+1][j] =  uniSetCorrectGenotype[cellID*2+1][snpUniID];
					k++;
				}
			}
			outputLgenFile(trueLgen, haplotypes);
		}
		File recodedtrueLgen = new File(workDir, ParamParser.outPrefix + "Set" + dataSetID + "_noErr.lgen");
		if(!recodedtrueLgen.exists())
			outputLgenFile(recodedtrueLgen, recode(haplotypes, tag20MPosTo2MIndex, recodeMap));
		return haplotypes;
	}

	public void readAlleleSwitchedSNPs(File switchedSNPFile, Map<Integer, Integer> tagPosToUniTagListIndex, boolean append) {
		if(!append || snpWithAlleleSwitched==null)
			snpWithAlleleSwitched = new HashSet<Integer>();
		FileHandler handleErrFile = new FileHandler(switchedSNPFile);
		handleErrFile.openForRead();
		String nextLine = null;
		while((nextLine = handleErrFile.nextLine()) != null)
			snpWithAlleleSwitched.add(tagPosToUniTagListIndex.get(Integer.valueOf(nextLine)));
		handleErrFile.closeReader();
	}

	public void readGenoShuffledSNPs(File shuffledFile, Map<Integer, Integer> tagPosToCurTagListIndex) {
		snpWithGenoShuffled = new HashSet<Integer>();
		FileHandler handleErrFile = new FileHandler(shuffledFile);
		handleErrFile.openForRead();
		String nextLine = null;
		while((nextLine = handleErrFile.nextLine()) != null)
			snpWithGenoShuffled.add(tagPosToCurTagListIndex.get(Integer.valueOf(nextLine)));
		handleErrFile.closeReader();	
	}

	public void readMissingErrs(File missingGenoFile) {
		missingGenoList = new HashMap<Integer, int[]>();
		FileHandler handleErrFile = new FileHandler(missingGenoFile);
		handleErrFile.openForRead();
		handleErrFile.nextLine();//header
		String nextLine = null;	
		while((nextLine = handleErrFile.nextLine()) != null){
			String[] s = nextLine.split("\t");
			int[] missingGeno = {Integer.valueOf(s[1]), Integer.valueOf(s[2])};
			missingGenoList.put(Integer.valueOf(s[0]), missingGeno);
		}
		handleErrFile.closeReader();
	}

	public void readGenoErrs(File wrongGenoFile) {
		wrongGenoList = new HashMap<Integer, int[]>();
		FileHandler handleErrFile = new FileHandler(wrongGenoFile);
		handleErrFile.openForRead();
		handleErrFile.nextLine();//header
		String nextLine = null;
		while((nextLine = handleErrFile.nextLine()) != null){
			String[] s = nextLine.split("\t");
			int[] wrongGeno = {Integer.valueOf(s[1]), Integer.valueOf(s[2]), Integer.valueOf(s[3]), Integer.valueOf(s[4]), Integer.valueOf(s[5]), Integer.valueOf(s[6])};
			wrongGenoList.put(Integer.valueOf(s[0]), wrongGeno);
		}
		handleErrFile.closeReader();
	}
	
	public void recordMissingErrs(File errorFile) {
	   	if(errorFile.exists()) errorFile.delete();
    	FileHandler handleErrFiles = new FileHandler(errorFile);
		handleErrFiles.openForWrite(true);
		handleErrFiles.writeString("Missing_Genotype_ID\tSubjectID\tMarker_Position");
		for(Entry<Integer, int[]> e : missingGenoList.entrySet()){
			int[] record = e.getValue();
			handleErrFiles.writeString(Integer.toString(e.getKey()) + "\t" + Integer.toString(record[0]) + "\t" + Integer.toString(record[1]));
		}
    	handleErrFiles.closeWriter();	
	}

	public void recordGenoErrs(File errorFile) {
    	//if(errorFile.exists()) errorFile.delete();
    	FileHandler handleErrFiles = new FileHandler(errorFile);
		handleErrFiles.openForWrite(true);
		handleErrFiles.writeString("Wrong_Genotype_ID\tSubjectID\tMarker_Position\tAlleles before change\tAlleles after change");
		for(Entry<Integer, int[]> e : wrongGenoList.entrySet()){
			int[] record = e.getValue();
			handleErrFiles.writeString(Integer.toString(e.getKey()) + "\t" + Integer.toString(record[0]) + "\t" + Integer.toString(record[1]) + "\t" 
					+ Integer.toString(record[2]) + "\t" + Integer.toString(record[3]) + "\t" + Integer.toString(record[4]) + "\t" + Integer.toString(record[5]));
		}	
    	handleErrFiles.closeWriter();
	}

	public void recordAlleleSwitchedSNPs(File errorFile) {
		//if(errorFile.exists()) errorFile.delete();
		FileHandler handleErrFiles = new FileHandler(errorFile);
		handleErrFiles.openForWrite(true);	
		for(Integer tagID : snpWithAlleleSwitched)
			handleErrFiles.writeString(Integer.toString(tagList.get(tagID)));
		handleErrFiles.closeWriter();
	}
	
	public Set<Integer> recordGenoShuffledSNPs(File errorFile) {
		Set<Integer> snp20Mpos = new HashSet<Integer>(snpWithGenoShuffled.size());
		if(errorFile.exists()) errorFile.delete();
		FileHandler handleErrFiles = new FileHandler(errorFile);
		handleErrFiles.openForWrite(true);
		for(Integer tagID : snpWithGenoShuffled){
			int snp20MPos = tagList.get(tagID);
			handleErrFiles.writeString(Integer.toString(snp20MPos));
			snp20Mpos.add(snp20MPos);
		}
		handleErrFiles.closeWriter();
		return snp20Mpos;
	}

	public void assignErr(int errType, int sizeOfPool, Set<Integer> itemToExclude, int numOfItemsToSelect, Map<Integer, Integer> tag20mPosTo2MIndex, List<Map<Character, Character>> recodeMap) {
		Set<Integer> selectedItems = new HashSet<Integer>();
		Random randGenerator = new Random();
		while(selectedItems.size() < numOfItemsToSelect){
			int i = randGenerator.nextInt(sizeOfPool);
			if(itemToExclude!=null && itemToExclude.contains(i)) continue;
			if(ParamParser.imposeErrOnATCGMarkerOnly && !isATCGMarker(recodeMap.get(tag20mPosTo2MIndex.get(tagList.get(i))))) continue;
			if(!ParamParser.imposeErrOnATCGMarkerOnly && isATCGMarker(recodeMap.get(tag20mPosTo2MIndex.get(tagList.get(i))))) continue;
			selectedItems.add(i);
		}
		List<Integer> sortedItems = new ArrayList<Integer>(selectedItems);
		selectedItems.clear();
		selectedItems = null;
		Collections.sort(sortedItems);
		switch(errType){
			case 0:	snpWithAlleleSwitched = new HashSet<Integer>(sortedItems);
					break;
			case 1:	snpWithGenoShuffled = new HashSet<Integer>(sortedItems);
					break;
			case 2: missingGenoList = new HashMap<Integer, int[]>(numOfItemsToSelect);
					for(int item : sortedItems)
						missingGenoList.put(item, null);
					break;
			case 3: wrongGenoList = new HashMap<Integer, int[]>(numOfItemsToSelect);
					for(int item : sortedItems)
						wrongGenoList.put(item, null);
					break;
		}
	}
		
	private boolean isATCGMarker(Map<Character, Character> map) {
		Set<Character> alleles = new HashSet<Character>(2);
		for(Character alleleState : map.values())
			alleles.add(alleleState);
		if(alleles.contains('A') && alleles.contains('T')) return true;
		if(alleles.contains('C') && alleles.contains('G')) return true;
		return false;
	}

}
