package seqSIM;

import java.io.File;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ErrController {
    public Map<Integer, int[]> allMissingGeno;//<subsetID, genotypeMatrixCellID>
    public Map<Integer, int[]> allWrongGeno;//<subsetID, genotypeMatrixCellID>
    public Set<Integer> allShuffledSNPs;
    
    public ErrController(){
    	if(ParamParser.missingRate > 0) allMissingGeno = new HashMap<Integer, int[]>();
    	if(ParamParser.errorRate > 0) allWrongGeno = new HashMap<Integer, int[]>();
    	if(ParamParser.randMarkerRate > 0) allShuffledSNPs = new HashSet<Integer>();
    }

	public void imposeError(String workDir, DataSet ds, Map<Integer, Integer> tagPosToCurrTagListIndex) {
		Integer SetID = (ParamParser.createSingleSet)? 1:ds.dataSetID;
		int tagCt = ds.getTagCt();
		if(ParamParser.switchRate>0){
			File switchedSNPFile = new File(workDir, ParamParser.outPrefix + "SwitchedSNPs_Set" + SetID.toString() + ".txt");
			if(switchedSNPFile.exists()) /*switchedSNPFile.delete();*/
				ds.readAlleleSwitchedSNPs(switchedSNPFile, tagPosToCurrTagListIndex, false);
			else{			
				ds.assignErr(0, tagCt, null, (int)(tagCt*ParamParser.switchRate));
				ds.recordAlleleSwitchedSNPs(switchedSNPFile);
			}
		}
		if(ParamParser.randMarkerRate > 0){
			File shuffledFile = new File(workDir, ParamParser.outPrefix + "ShuffledSNPs_Set" + SetID.toString() + ".txt");			
			if(shuffledFile.exists())/*shuffledFile.delete();*/
				ds.readGenoShuffledSNPs(shuffledFile, tagPosToCurrTagListIndex);
			else{
				ds.assignErr(1, tagCt, ds.snpWithAlleleSwitched, (int)(tagCt*ParamParser.randMarkerRate));
				allShuffledSNPs.addAll(ds.recordGenoShuffledSNPs(shuffledFile));	
			}
		}
		if(ParamParser.missingRate > 0 || ParamParser.errorRate > 0){
			int genoCt = ds.getDataSetSize()*ds.getTagCt();
			if(ParamParser.missingRate > 0){
				File missingGenoFile = new File(workDir, ParamParser.outPrefix + "MissingGeno_Set" + SetID.toString() + ".txt");
				if(missingGenoFile.exists())/*missingGenoFile.delete();*/
					ds.readMissingErrs(missingGenoFile);
				else
					ds.assignErr(2, genoCt, null,(int)(genoCt*ParamParser.missingRate));
			}
			if(ParamParser.errorRate > 0){
				File wrongGenoFile = new File(workDir, ParamParser.outPrefix + "WrongGeno_Set" + SetID.toString() + ".txt");
				if(wrongGenoFile.exists()) /*wrongGenoFile.delete();*/
					ds.readGenoErrs(wrongGenoFile);
				else
					ds.assignErr(3, genoCt, ds.missingGenoList.keySet(), (int)(genoCt* ParamParser.errorRate));
			}
		}	
	}

/*	public void imposeTestingErrOnSingleSet(String workDir, DataSet curDataSet, Map<Integer, Integer> tagPosToUniTagListIndex) {
		//TODO: so far it only impose flipping errors by first finding the union of flipped snps in all subsets and then flipping these snps in the uniset.
		File switchedSNPInCurDataSet = new File(workDir, ParamParser.outPrefix + "SwitchedSNPs_Set0.txt");
		if(switchedSNPInCurDataSet.exists()) switchedSNPInCurDataSet.delete();
			curDataSet.readAlleleSwitchedSNPs(switchedSNPInCurDataSet, tagPosToUniTagListIndex, false);
		else{
			curDataSet.snpWithAlleleSwitched = null;
			for(int setID = 1; setID <= ParamParser.num_subsets; setID++){
				File switchedSNPInSubSet = new File(workDir, ParamParser.outPrefix + "SwitchedSNPs_Set" + setID + ".txt");
				if(switchedSNPInSubSet.exists())
					curDataSet.readAlleleSwitchedSNPs(switchedSNPInSubSet, tagPosToUniTagListIndex, true);
			}
			curDataSet.recordAlleleSwitchedSNPs(switchedSNPInCurDataSet);
		}
	}*/

	public void imposeError(String workDir, DataSet ds, Map<Integer, Integer> tagPosToCurrTagListIndex, Map<Integer, Integer> tag20mPosTo2MIndex,
			List<Map<Character, Character>> recodeMap) {
		Integer SetID = (ParamParser.createSingleSet)? 1:ds.dataSetID;
		int tagCt = ds.getTagCt();
		if(ParamParser.switchRate > 0){
			File switchedSNPFile = new File(workDir, ParamParser.outPrefix + "SwitchedSNPs_Set" + SetID.toString() + ".txt");
			if(switchedSNPFile.exists())/*switchedSNPFile.delete();*/
				ds.readAlleleSwitchedSNPs(switchedSNPFile, tagPosToCurrTagListIndex, false);
			else{
				int numOfFlippedMarker = (ParamParser.singleFlip)?1:(int)(tagCt*ParamParser.switchRate);
				if(ParamParser.imposeErrBasedOAlleleStates)
					ds.assignErr(0, tagCt, null, numOfFlippedMarker,tag20mPosTo2MIndex,recodeMap);
				else
					ds.assignErr(0, tagCt, null, numOfFlippedMarker);
				ds.recordAlleleSwitchedSNPs(switchedSNPFile);
			}
		}
		if(ParamParser.randMarkerRate > 0){
			File shuffledFile = new File(workDir, ParamParser.outPrefix + "ShuffledSNPs_Set" + SetID.toString() + ".txt");
			if (shuffledFile.exists()) {
				/*shuffledFile.delete();*/
				ds.readGenoShuffledSNPs(shuffledFile, tagPosToCurrTagListIndex);
				for (int curTagListID : ds.snpWithGenoShuffled)
					allShuffledSNPs.add(ds.tagList.get(curTagListID));
			} else {
				if (!ParamParser.completeRandShuffle && ParamParser.startFromHomozygote && ParamParser.focusOnMinorHomo){
				}else{
					if (ParamParser.imposeErrBasedOAlleleStates)
						ds.assignErr(1, tagCt, ds.snpWithAlleleSwitched, (int) (tagCt * ParamParser.randMarkerRate),
								tag20mPosTo2MIndex, recodeMap);
					else
						ds.assignErr(1, tagCt, ds.snpWithAlleleSwitched, (int) (tagCt * ParamParser.randMarkerRate));
					allShuffledSNPs.addAll(ds.recordGenoShuffledSNPs(shuffledFile));
				}
			} 
		}
		if(ParamParser.missingRate > 0 || ParamParser.errorRate > 0){
			int genoCt = ds.getDataSetSize()*ds.getTagCt();
			if(ParamParser.missingRate > 0){
				File missingGenoFile = new File(workDir, ParamParser.outPrefix + "MissingGeno_Set" + SetID.toString() + ".txt");
				if(missingGenoFile.exists()) /*missingGenoFile.delete();*/
					ds.readMissingErrs(missingGenoFile);
				else
					ds.assignErr(2, genoCt, null,(int)(genoCt*ParamParser.missingRate));
			}
			if(ParamParser.errorRate > 0){
				File wrongGenoFile = new File(workDir, ParamParser.outPrefix + "WrongGeno_Set" + SetID.toString() + ".txt");
				if(wrongGenoFile.exists()) /*wrongGenoFile.delete();*/
					ds.readGenoErrs(wrongGenoFile);
				else
					ds.assignErr(3, genoCt, ds.missingGenoList.keySet(), (int)(genoCt* ParamParser.errorRate));
			}
		}	
		
	}				
}
