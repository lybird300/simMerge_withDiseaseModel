package seqSIM;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Scanner;
import java.util.Set;

import haploview.tagger.TaggerController;

/**
 * A tool for creating a human disease GWAS data set using DNA sequences generated from a coalescent simulation
 * !!!!!!!!!!!!!!!!!!!!Alert: this version of code can only work on a largemem node, as it directly inputs coJava hap file
 * @author Yuan Lin
 */
public class SimMergeMain {
	private static ParamParser params;
    //Mapping a SNP to the frequency of its mutation allele; entry format: <index of a SNP regarding the all-SNP array (stored in out.snp-1.gz), its MAF regarding a specific panel> 
    private static Map<Integer, Double> mutAlleleFreq = new HashMap<Integer, Double>();
    
    //Mapping a SNP to chromosomes with mutation on that SNP; entry format: <index of a SNP regarding the all-SNP array (stored in out.snp-1.gz), chromID>
    private static Map<Integer, Set<Integer>> snpToChrom = new HashMap<Integer, Set<Integer>>();
    
	
	public static void main(String[] args) throws IOException, InterruptedException {
		if(args.length == 0){
            System.out.println("\nPlease indicate name of the .properties file that provides default values for program parameters.\n");
            System.exit(0);
		}
        else{
	       	// Obtain and assign parameter values
	       	setParams(new ParamParser(args[0]));
	       	
	       	//create the 1st level directory that holds all data output	
	    	String dir_1st_name = ParamParser.dataOutputDir /*+ File.separator + "simGWASNewData"*/;
	    	File dir_1st = new File(dir_1st_name);
	       	if(!dir_1st.exists()) dir_1st.mkdirs();
		
	       	//check selected casual variants (at the first level directory)
			File cvFile = new File(dir_1st_name, "selectedCausalVariants.txt");
			FileHandler handleCausalVar = null;
			//If the file does not exist, it means causal variants have not been selected yet. Go ahead selecting them and write relevant 
	      	//information into the file. The minorAlleleFreq map and the snpToChrom map will be created during the process. 	     
			if(!cvFile.exists()){
				handleCausalVar = new FileHandler(cvFile);
				handleCausalVar.openForWrite(true);
				handleCausalVar.writeString("CV Index\tMAF\tChromosomes with this CV");
				//Random r = new Random();
				double cvMafLowerBound = ParamParser.cvMafSets[0];
				double cvMafUpperBound = ParamParser.cvMafSets[1];
				//double minorAlleleFreq = cvMafLowerBound+r.nextDouble()*(cvMafUpperBound-cvMafLowerBound);
				selectCausalVariants(ParamParser.segmentStart, ParamParser.segmentEnd, cvMafLowerBound, cvMafUpperBound, ParamParser.numOfCVPerMAFRange);
				for(Entry<Integer, Double> var : mutAlleleFreq.entrySet()){
					int cv = var.getKey();
				    List<Integer> chromsWithMut = new ArrayList<Integer>(snpToChrom.get(cv));
				    String chromLine = "";
				    for(Integer chrom : chromsWithMut)
				      	chromLine += "\t" + chrom;
				    double cvMaf = (var.getValue()>0.5)?(1-var.getValue()):var.getValue();
			      	handleCausalVar.writeString(String.valueOf(cv)+ "\t" +  String.valueOf(cvMaf) + chromLine);
				}
	      		handleCausalVar.closeWriter();
	      		System.exit(0);
			}
			
	      	// If the file already exists, it means causal variants have been selected and we are working in a recovery mode.
			// Then read their info from the file to populate the minorAlleleFreq map and the snpToChrom map. 
	      	handleCausalVar = new FileHandler(cvFile);
	      	handleCausalVar.openForRead();
	      	handleCausalVar.nextLine();//header
			String cvRecord = null;
			while((cvRecord = handleCausalVar.nextLine())!=null){
				Scanner cvLine = new Scanner(cvRecord);
				cvLine.useDelimiter("\t");
				int aSNP = Integer.parseInt(cvLine.next());
				mutAlleleFreq.put(aSNP, Double.parseDouble(cvLine.next()));
				while(cvLine.hasNext()){
					Set<Integer> chromsWithSpecificSNPMut = snpToChrom.get(aSNP);			
					if (chromsWithSpecificSNPMut == null) {
						chromsWithSpecificSNPMut = new HashSet<Integer>();
						snpToChrom.put(aSNP, chromsWithSpecificSNPMut);
					}
					chromsWithSpecificSNPMut.add(Integer.parseInt(cvLine.next()));
				}
				cvLine.close();
			}
			handleCausalVar.closeReader();
	      	
			int cv = Integer.valueOf(args[1]);
			Set<Integer> chromsWithSNP = snpToChrom.get(cv);
			double cvMaf = mutAlleleFreq.get(cv);			
			mutAlleleFreq = null;
			snpToChrom = null;
			/*int caseIndex = Integer.valueOf(args[2]);
			int ctrlIndex = Integer.valueOf(args[3]);*/
			int grrIndex = Integer.valueOf(args[2]);
			int bdrIndex = Integer.valueOf(args[3]);
			int repIndex = Integer.valueOf(args[4]);
			List<Integer> chromsWithCV = new ArrayList<Integer>(chromsWithSNP);
			List<Integer> chromsWithoutCV = new ArrayList<Integer>();
	  		for(int k=0; k<ParamParser.sampleSize; k++){
	  			if(!chromsWithCV.contains((Integer)k))
	  				chromsWithoutCV.add((Integer)k);
	  		}
	  		
	      	// At the 2nd level directory, create a folder for the current SNP and name the folder as:
			// SNP_pos<its position>_freqLo<lower frequency bound>_freqHi<higher frequency bound>
			File dir_2nd = Util.returnSNPDir(cv, cvMaf);
			if(!dir_2nd.exists()) dir_2nd.mkdirs();
			/* The following code outputs all data related to the first replicate of a case-control panel.
			 * The chromosomes selected for different replicates of the case-control panel is stored in a file called "ChromsUsed***.txt"
			 * inside the folder for that specific panel and will be used later to figure out the SNP markers for that specific panel. */	
			CaseControlSampler sampler = new CaseControlSampler(dir_2nd, grrIndex, bdrIndex, repIndex);
			String workDir = sampler.getLocation();
			List<DataSet> dataSetList = sampler.prepPanel(cvMaf, chromsWithCV, chromsWithoutCV);
			System.out.println("The .fam files are ready.");
			chromsWithCV.clear();
			chromsWithCV = null;
						
			Map<Integer, Integer> allsnpPosToIndex = new HashMap<Integer, Integer>();			
			File refMapFile = new File(dir_2nd, ParamParser.outPrefix + "reference.map");
			boolean refMapExist = refMapFile.exists();
			FileHandler writeRefMapFile = null;
	  		if(!refMapExist){
	  			writeRefMapFile = new FileHandler(refMapFile);
		  		writeRefMapFile.openForWrite(true);
	  		}
	  		List<Map<Character, Character>> recodeMap = new ArrayList<Map<Character, Character>>();
	  		File alleleRecodeMap = new File(dir_2nd, ParamParser.outPrefix + "alleleCoding.gz");
	  		boolean recodeMapExist = alleleRecodeMap.exists();
	  		FileHandler readSNPPos = null;
	  		FileHandler writeRecodeMap = null;
	  		if(recodeMapExist)
	  			readSNPPos = new FileHandler(alleleRecodeMap);
	  		else{
	  			writeRecodeMap = new FileHandler(alleleRecodeMap);
	  			writeRecodeMap.openForWrite(true);
	  			readSNPPos = new FileHandler(ParamParser.dataInputDir, ParamParser.posInput);
	  		}
			readSNPPos.openForRead();
			String nextLine = null;	
			int snp2MIndx = 0;
			while((nextLine = readSNPPos.nextLine()) != null){
				String[] s = nextLine.split("\\s+");
				int snp20MPos = Integer.valueOf(s[0]);
				allsnpPosToIndex.put(snp20MPos, snp2MIndx);
				if(!refMapExist) writeRefMapFile.writeString(ParamParser.founderChrom + "\t" + "rs" + snp20MPos + "\t0\t" + snp20MPos);
				Character[] allele4Type = null;
				if(!recodeMapExist){
					allele4Type = Util.randomSelectAlleleWithoutReplacement();
					writeRecodeMap.writeString(snp20MPos + "\t" + allele4Type[0] + "\t" + allele4Type[1]);
				}
				else{
					allele4Type = new Character[2];
					allele4Type[0] = s[1].charAt(0);
					allele4Type[1] = s[2].charAt(0);
				}
				Map<Character, Character> allele2TypeTo4Type = new HashMap<Character, Character>(2);
				allele2TypeTo4Type.put('1', allele4Type[0]);
				allele2TypeTo4Type.put('2', allele4Type[1]);
				recodeMap.add(allele2TypeTo4Type);
				snp2MIndx++;
			}
			readSNPPos.closeReader();
			if(!refMapExist) writeRefMapFile.closeWriter();
			if(!recodeMapExist) writeRecodeMap.closeWriter();
			
			createRefHaps(dir_2nd, chromsWithoutCV, sampler.getUsedChroms(), ParamParser.founderChrom, recodeMap);
			System.out.println("References are ready.");
			
			Map<Integer, Integer> tag20MPosTo2MIndex = new HashMap<Integer, Integer>();
			TaggerController tagFinder = new TaggerController();
		    Map<Integer, List<Integer>> tagSNPs = tagFinder.prepTagSNPs(workDir, cv, chromsWithoutCV, allsnpPosToIndex, tag20MPosTo2MIndex);
		    if(ParamParser.createSingleSet){
		    	dataSetList.get(0).tagList = tagFinder.getUniTagList();
		    }
		    else{
		    	for(int subsetID = 1; subsetID <= ParamParser.num_subsets; subsetID++)
		    		dataSetList.get(subsetID).tagList = tagSNPs.get(subsetID);
		    	dataSetList.get(0).tagList = tagFinder.getUniTagList();
		    }
		    chromsWithoutCV.clear();
		    chromsWithoutCV = null;
		    allsnpPosToIndex.clear();
		    allsnpPosToIndex = null;
   			System.out.println("The .map files are ready.");
		   
		    Map<Integer, Map<Integer,Integer>> tagPosToTagListIndex = new HashMap<Integer, Map<Integer,Integer>>(ParamParser.num_subsets+1);
		    for(int setID = 0; setID <= ParamParser.num_subsets; setID++){
		    	tagPosToTagListIndex.put(setID, new LinkedHashMap<Integer, Integer>());
		    	Collections.sort(dataSetList.get(setID).tagList);
		    	for(int tID = 0; tID < dataSetList.get(setID).tagList.size(); tID++)
		    		tagPosToTagListIndex.get(setID).put(dataSetList.get(setID).tagList.get(tID), tID);
		    }
		    
		    ErrController ec = new ErrController();
		    if(ParamParser.createSingleSet){
		    	DataSet ds = dataSetList.get(0);
		    	if(ParamParser.imposeErrBasedOAlleleStates)
		    		ec.imposeError(workDir, ds, tagPosToTagListIndex.get(0),tag20MPosTo2MIndex,recodeMap);
		    	else
		    		ec.imposeError(workDir, ds, tagPosToTagListIndex.get(0));
		    }
		    else{
			    for(int setID = 1; setID < dataSetList.size(); setID++){			    	
			    	DataSet ds = dataSetList.get(setID);
			    	if(ParamParser.imposeErrBasedOAlleleStates)
			    		ec.imposeError(workDir, ds, tagPosToTagListIndex.get(setID),tag20MPosTo2MIndex,recodeMap);
			    	else
			    		ec.imposeError(workDir, ds, tagPosToTagListIndex.get(setID));
			    }
		    }
		    /*if(ParamParser.imposeFlipErrOnUniSet)
		    	ec.imposeTestingErrOnSingleSet(workDir, dataSetList.get(0), tagPosToUniTagListIndex);*/
		    Integer[] usedChromsArray = sampler.getUsedChroms().toArray(new Integer[sampler.getUsedChroms().size()]);
		    int[][] subIDtoHapID = sampler.getsubHapID();
		    byte[][] uniSetCorrect12Haplotypes = dataSetList.get(0).getUniSetCorrect12Haplotypes(workDir, tagPosToTagListIndex.get(0), 
		    		usedChromsArray, tag20MPosTo2MIndex, recodeMap, subIDtoHapID);
		    subIDtoHapID = null;
		    usedChromsArray = null;
		    if(!ParamParser.createSingleSet){
			    for(int setID = 1; setID < dataSetList.size(); setID++){
			    	DataSet ds = dataSetList.get(setID);
			    	byte[][] subsetCorrect12Haplotypes = ds.getSubSetCorrect12Haplotypes(workDir, tagPosToTagListIndex.get(0), 
			    				uniSetCorrect12Haplotypes, tag20MPosTo2MIndex, recodeMap);
			    	ds.generateSubSetErrHaplotypes(workDir, subsetCorrect12Haplotypes, tag20MPosTo2MIndex, recodeMap, ec);
			    	if(ParamParser.errorRate > 0){
			    		File wrongGenoFile = new File(workDir, ParamParser.outPrefix + "WrongGeno_Set" + ds.dataSetID + ".txt");
			    		if(!wrongGenoFile.exists()) ds.recordGenoErrs(wrongGenoFile);
			    	}
			    	if(ParamParser.missingRate > 0){
			    		File missingGenoFile = new File(workDir, ParamParser.outPrefix + "MissingGeno_Set" + ds.dataSetID + ".txt");
			    		if(!missingGenoFile.exists()) ds.recordMissingErrs(missingGenoFile);
			    	}
					if(ds.missingGenoList!= null) ec.allMissingGeno.putAll(ds.missingGenoList);
					if(ds.wrongGenoList != null) ec.allWrongGeno.putAll(ds.wrongGenoList);
			    }
		    }
		    if(ParamParser.missingRate > 0 || ParamParser.errorRate > 0 || ParamParser.randMarkerRate > 0 || ParamParser.imposeFlipErrOnUniSet || ParamParser.createSingleSet)
		    	dataSetList.get(0).generateUniSetErrHaplotypes(workDir, tagPosToTagListIndex.get(0), uniSetCorrect12Haplotypes, ec, tag20MPosTo2MIndex, recodeMap);
		    tag20MPosTo2MIndex.clear();
		    tag20MPosTo2MIndex = null;
		    recodeMap.clear();
		    recodeMap = null;
		    System.out.println("The .lgen files are ready.");
		    //System.out.println("Note: Errors in the uniset is a combination of errors in all subsets, so there is no corresponding error files for Set0.");	    
        }
	}

	public static void createRefHaps(File workDir, List<Integer> chromsWithoutCV, Set<Integer> usedChroms, int chrID, List<Map<Character, Character>> recodeMap){
		File refFile = new File(workDir, ParamParser.outPrefix + "reference_chr" + chrID + ".gz");
		if(refFile.exists()) return;
		File oldRefFile = new File(workDir, ParamParser.outPrefix + "reference_chr" + chrID + "_allele12.gz");
		char[][] hapList = new char[ParamParser.num_refSub*2][];
		if(!oldRefFile.exists()){
			chromsWithoutCV.removeAll(usedChroms);
			System.out.println(chromsWithoutCV.size() + " hap samples remain.");
	  		Collections.shuffle(chromsWithoutCV);
	  		Integer[] refHaps = new Integer[ParamParser.num_refSub*2];
			for(int j = 0; j < refHaps.length; j++)
				refHaps[j] = chromsWithoutCV.get(j);
			Arrays.sort(refHaps);
			Collection<String> hapStrings = Util.getHapsfromDistFiles(refHaps, null).values();
			int sampleHapCt = 0;
			for(String hap : hapStrings){
				hapList[sampleHapCt++] = hap.toCharArray();
			}
		}
		else{
			FileHandler readOldRefFile = new FileHandler(oldRefFile);
			readOldRefFile.openForRead();
			String nextLine = null;
			int sampleHapCt = 0;
			while((nextLine=readOldRefFile.nextLine())!=null){
				String[] s = nextLine.split("\\s+");
				hapList[sampleHapCt++] = s[2].toCharArray();
			}
			readOldRefFile.closeReader();
		}
		for(int snpCt = 0; snpCt < recodeMap.size(); snpCt++){
			Map<Character, Character> snpMap = recodeMap.get(snpCt);
			for(int hapCt = 0; hapCt < hapList.length; hapCt++){
				Character oldAllele = hapList[hapCt][snpCt];
				hapList[hapCt][snpCt] = snpMap.get(oldAllele);
			}
		}
		FileHandler writeRefFile = new FileHandler(refFile);
  		writeRefFile.openForWrite(true);
  		int sampleHapCt = 0;
		for(;sampleHapCt < hapList.length; sampleHapCt++){			
			String ref = "Sub" + sampleHapCt/2;
			if(sampleHapCt < 2)
				ref += (sampleHapCt==0)?"\tHAPLO1\t":"\tHAPLO2\t";
			else
				ref += (sampleHapCt%2==0)?"\tHAPLO1\t":"\tHAPLO2\t";
			ref += String.valueOf(hapList[sampleHapCt]);
			writeRefFile.writeString(ref);
			ref = null;
		}	
  		writeRefFile.closeWriter();
	}
	
/*	private static String recode(String oldSeq, List<Map<Byte, Character>> recodeMap) {
		char[] seqArray = oldSeq.toCharArray();
		String newSeq = new String("");
		for(int i = 0; i < seqArray.length; i++)
			newSeq += recodeMap.get(i).get(Byte.valueOf(String.valueOf(seqArray[i])));
		return newSeq;
	}*/

	private static void selectCausalVariants(int pos_start, int pos_end, double cvMafLowerBound, double cvMafUpperBound, int maxNum) {	
		Map<String, ArrayList<Integer>> chromRecord = Util.getOriginDataLocation(ParamParser.mutInput, null);
		int ctReadPosFile = 0;
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){
			ctReadPosFile++;
			if(ctReadPosFile==ParamParser.inputFileNum)
				System.out.println("Last pos file.\n");
			FileHandler getMutationData = new FileHandler(ParamParser.dataInputDir, entry.getKey());
			getMutationData.openForRead();
			getMutationData.nextLine();//first line is the header			
			String nextLine = null;
			Integer mutInx;
			String[] s = null;
			while((nextLine = getMutationData.nextLine()) != null){
				if(nextLine=="" || nextLine=="\n") continue;
				s = nextLine.split("\\s+");
				int chromID = Integer.parseInt(s[0]);
				for(int j=1; j<s.length; j++){
					mutInx = Integer.parseInt(s[j]);
					boolean exceedFreqThresh = false;
					
					//update the minorAlleleFreq map
					if(mutAlleleFreq.get(mutInx)!=null){
						double currentValue = mutAlleleFreq.get(mutInx).doubleValue() + 1/(double)ParamParser.sampleSize;
						mutAlleleFreq.put(mutInx, currentValue);						
						if(currentValue > 0.5 && currentValue > 1-cvMafUpperBound)
							exceedFreqThresh = true;
					}
					else
						mutAlleleFreq.put(mutInx, 1/(double)ParamParser.sampleSize);
					
					if(!exceedFreqThresh){//update the snpToChrom map only when exceedFreqThresh == false
						Set<Integer> chromsWithSpecificSNPMut = snpToChrom.get(mutInx);			
						if (chromsWithSpecificSNPMut == null) {
							//initiate the set of chromosomes with mutation on a specific snp, if it does not exist
							//The Set data structure is used to avoid redundancy
							chromsWithSpecificSNPMut = new HashSet<Integer>();//
							snpToChrom.put(mutInx, chromsWithSpecificSNPMut);
						}
						chromsWithSpecificSNPMut.add(chromID);		
					}
					else if(snpToChrom.get(mutInx) != null)
						snpToChrom.remove(mutInx);
				}
				
			}
			getMutationData.closeReader();
		}		
		/*System.out.println("In this data set, MAF ranges from " + Collections.min(minorAlleleFreq.values())
				+ " to " + Collections.max(minorAlleleFreq.values()) + ".\n");*/
		List<Integer> allSNPs = new ArrayList<Integer>(mutAlleleFreq.keySet());
	    Collections.shuffle(allSNPs);
	    //System.out.println("no\n");
		Integer[] snpPos = new Integer[ParamParser.num_mut];
		//Read the absolute positions of all SNPs from the "out.snp-1.gz" file
		FileHandler readSNPPos = new FileHandler(ParamParser.dataInputDir, ParamParser.posInput);
		readSNPPos.openForRead();
		readSNPPos.nextLine();//header
		String nextLine = null;
		int j = 0;
		while((nextLine = readSNPPos.nextLine()) != null)
			snpPos[j++] = Integer.valueOf(nextLine);
		readSNPPos.closeReader();
	    
	    Map<Integer, Double> selectedSNPs_maf = new HashMap<Integer, Double>();
	    Map<Integer, Set<Integer>> selectedSNPs_chroms = new HashMap<Integer, Set<Integer>>();
		for (int i = 0; i < allSNPs.size(); i++){
	      	int snpIndex = allSNPs.get(i);
	      	double minorAlleleFreq = mutAlleleFreq.get(snpIndex);
	      	if(minorAlleleFreq > 0.5) minorAlleleFreq = 1-minorAlleleFreq;
	      	int realPos = snpPos[snpIndex];
			if (realPos < pos_start || realPos > pos_end || minorAlleleFreq < cvMafLowerBound || minorAlleleFreq > cvMafUpperBound)
				continue;
			Set<Integer> c = snpToChrom.get(snpIndex);
			if(c == null)
				continue;
			else{
	      		selectedSNPs_maf.put(realPos, minorAlleleFreq);
	      		selectedSNPs_chroms.put(realPos, c);
			}
			if(selectedSNPs_maf.size() >= maxNum) break;
		}
		mutAlleleFreq = selectedSNPs_maf;
		snpToChrom = selectedSNPs_chroms;
	}

	public static ParamParser getParams() {
		return params;
	}

	public static void setParams(ParamParser params) {
		SimMergeMain.params = params;
	}
	
	


}
