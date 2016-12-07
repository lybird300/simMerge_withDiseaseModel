package seqSIM;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.FieldPosition;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Locale;
import java.util.Map;
import java.util.Random;
import java.util.Map.Entry;

import seqSIM.FileHandler;
import seqSIM.ParamParser;
import seqSIM.Subject;

public abstract class Util {

    public static String formatPValue(double pval){
        DecimalFormat df;
        //java truly sucks for simply restricting the number of sigfigs but still
        //using scientific notation when appropriate
        if (pval < 0.0001){
            df = new DecimalFormat("0.0000E0", new DecimalFormatSymbols(Locale.US));
        }else{
            df = new DecimalFormat("0.0000", new DecimalFormatSymbols(Locale.US));
        }
        String formattedNumber =  df.format(pval, new StringBuffer(), new FieldPosition(NumberFormat.INTEGER_FIELD)).toString();
        return formattedNumber;
    }

    public static double roundDouble (double d, int places){
        double factor = Math.pow(10,places);
        return Math.rint(d*factor)/factor;
    }
    
	//Get indices of specific marker SNPs on selected chromsomes.
	public static Map<Integer, byte[]> gatherMutMarkersfromDistFiles(Integer[] chroms, Integer[] markerIndex) {		
		Arrays.sort(chroms);
		Arrays.sort(markerIndex);
		Map<Integer, byte[]> chromsToMarkerSNPs = new HashMap<Integer, byte[]>();
		Map<String, ArrayList<Integer>> chromRecord = getOriginDataLocation(ParamParser.mutInput, chroms);	
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){
			FileHandler readMarkerData = new FileHandler(ParamParser.dataInputDir, entry.getKey());
			readMarkerData.openForRead();
			readMarkerData.nextLine();//first line is the header
			ArrayList<Integer> records = entry.getValue();
			int maxRecIndex = Collections.max(records);
			Integer[] recordsInThisFile = records.toArray(new Integer[records.size()]);
			Arrays.sort(recordsInThisFile);		
			String nextLine = null;
			for(int recordCt = 0; (nextLine = readMarkerData.nextLine()) != null && recordCt <= maxRecIndex; recordCt++){
				if(Arrays.binarySearch(recordsInThisFile, recordCt) < 0)
					continue;//this record is not what we want
				String[] s = nextLine.split("\\s+");					
				int chromID = Integer.parseInt(s[0]);//In out.pos-1_# files, s[0] is chromID
				if(Arrays.binarySearch(chroms, chromID) < 0)
					System.out.println("There is something wrong with you strategy in reading distributed files.\n");
				int[] mutations = new int[s.length-1];
				for(int j = 0; j < s.length-1; j++)
					mutations[j] = Integer.valueOf(s[j+1]);
				Arrays.sort(mutations);
				byte[] markers = new byte[markerIndex.length];
				for(int k = 0; k < markerIndex.length; k++){
					if(Arrays.binarySearch(mutations, markerIndex[k]) >= 0)
						markers[k] = (byte)1;
					else
						markers[k] = (byte)2;
				}
				chromsToMarkerSNPs.put(chromID, markers);				
			}				
			readMarkerData.closeReader();
		}
		return chromsToMarkerSNPs;
	}
	
	public static int[] gatherMutMarkersfromDistFiles(int chromID, Integer[] markerIndex) {		
		Arrays.sort(markerIndex);
		int[] markers = new int[markerIndex.length];
		int fileSize = ParamParser.sampleSize/ParamParser.inputFileNum;
		String fileName = ParamParser.mutInput + "_" + (chromID/fileSize) + ".gz";
		FileHandler readMarkerData = new FileHandler(ParamParser.dataInputDir, fileName);
		readMarkerData.openForRead();
		readMarkerData.nextLine();//first line is the header		
		String nextLine = null;
		while((nextLine = readMarkerData.nextLine()) != null){
			String[] s = nextLine.split("\\s+");					
			if(Integer.parseInt(s[0]) != chromID)//In out.pos-1_# files, s[0] is chromID
				continue;
			int[] mutations = new int[s.length-1];
			for(int j = 0; j < s.length-1; j++)
				mutations[j] = Integer.valueOf(s[j+1]);
			Arrays.sort(mutations);
			
			for(int k = 0; k < markerIndex.length; k++){
				if(Arrays.binarySearch(mutations, markerIndex[k]) >= 0)
					markers[k] = 1;
				else
					markers[k] = 2;
			}
			break;
		}
		readMarkerData.closeReader();
		return markers;
	}
	
	public static char[] getFullHap(int chromID) {		
		char[] markers = new char[ParamParser.num_mut];
		int fileSize = ParamParser.sampleSize/ParamParser.inputFileNum;
		String fileName = ParamParser.mutInput + "_" + (chromID/fileSize) + ".gz";
		FileHandler readMarkerData = new FileHandler(ParamParser.dataInputDir, fileName);
		readMarkerData.openForRead();
		readMarkerData.nextLine();//first line is the header		
		String nextLine = null;
		while((nextLine = readMarkerData.nextLine()) != null){
			String[] s = nextLine.split("\\s+");					
			if(Integer.parseInt(s[0]) != chromID)//In out.pos-1_# files, s[0] is chromID
				continue;
			int[] mutations = new int[s.length-1];
			for(int j = 0; j < s.length-1; j++)
				mutations[j] = Integer.valueOf(s[j+1]);
			Arrays.sort(mutations);
			
			for(int k = 0; k < markers.length; k++){
				if(Arrays.binarySearch(mutations, k) >= 0)
					markers[k] = '1';
				else
					markers[k] = '2';
			}
			break;
		}
		readMarkerData.closeReader();
		return markers;
	}
	
	public static Map<Integer, String> getHapsfromDistFiles(Integer[] chroms, Integer[] markerIndex) {		
		Arrays.sort(chroms);
		if(markerIndex!=null) Arrays.sort(markerIndex);
		Map<Integer, String> chromsToMarkerSNPs = new HashMap<Integer, String>();
		Map<String, ArrayList<Integer>> chromRecord = getOriginDataLocation(ParamParser.mutInput, chroms);	
		for(Entry<String, ArrayList<Integer>> entry: chromRecord.entrySet()){
			FileHandler readMarkerData = new FileHandler(ParamParser.dataInputDir, entry.getKey());
			readMarkerData.openForRead();
			readMarkerData.nextLine();//first line is the header
			ArrayList<Integer> records = entry.getValue();
			int maxRecIndex = Collections.max(records);
			Integer[] recordsInThisFile = records.toArray(new Integer[records.size()]);
			Arrays.sort(recordsInThisFile);		
			String nextLine = null;
			for(int recordCt = 0; (nextLine = readMarkerData.nextLine()) != null && recordCt <= maxRecIndex; recordCt++){
				if(Arrays.binarySearch(recordsInThisFile, recordCt) < 0)
					continue;//this record is not what we want
				String[] s = nextLine.split("\\s+");					
				int chromID = Integer.parseInt(s[0]);//In out.pos-1_# files, s[0] is chromID
				if(Arrays.binarySearch(chroms, chromID) < 0)
					System.out.println("There is something wrong with you strategy in reading distributed files.\n");
				int[] mutations = new int[s.length-1];
				for(int j = 0; j < s.length-1; j++)
					mutations[j] = Integer.valueOf(s[j+1]);
				Arrays.sort(mutations);
				int length, compareValue;
				if(markerIndex==null) length = ParamParser.num_mut;
				else length = markerIndex.length; 
				char[] markers = new char[length];
				for(int k = 0; k < length; k++){
					if(markerIndex==null) compareValue = k;
					else compareValue = markerIndex[k];
					if(Arrays.binarySearch(mutations, compareValue) >= 0)
						markers[k] = '1';
					else
						markers[k] = '2';
				}
				chromsToMarkerSNPs.put(chromID, new String(markers));				
			}				
			readMarkerData.closeReader();
		}
		return chromsToMarkerSNPs;
	}
    
	/**
	 * Returns a pseudo-random number between min and max, inclusive.
	 * The difference between min and max can be at most
	 * <code>Integer.MAX_VALUE - 1</code>.
	 *
	 * @param min Minimum value
	 * @param max Maximum value.  Must be greater than min.
	 * @return Integer between min and max, inclusive.
	 * @see java.util.Random#nextInt(int)
	 */
	public static int randInt(int min, int max) {
		Random randomGenerator = new Random();
	    // nextInt is normally exclusive of the top value, so add 1 to make it inclusive
	    return randomGenerator.nextInt((max - min) + 1) + min;
	}
	
	/**
	 * Find out where the information of each chromosome exists and store the results in a map.
	 * Each entry of this map is <name of the file, line in the file>
	 * @param filePrefix
	 * @param usedChroms
	 * @return
	 */
	public static Map<String, ArrayList<Integer>> getOriginDataLocation(String filePrefix, Integer[] usedChroms){
		Map<String, ArrayList<Integer>> fileNameToChromIDs = new HashMap<String, ArrayList<Integer>>();
		if(usedChroms==null){
			for(int j = 0; j< ParamParser.inputFileNum; j++)
				fileNameToChromIDs.put(filePrefix + "_" + j + ".gz", new ArrayList<Integer>());
		}
		else{
			for(int i = 0; i<usedChroms.length; i++){
				int fileSize = ParamParser.sampleSize/ParamParser.inputFileNum;
				String fileName = filePrefix + "_" + (usedChroms[i]/fileSize) + ".gz";
				if(fileNameToChromIDs.get(fileName)==null){
					fileNameToChromIDs.put(fileName, new ArrayList<Integer>());	
				}
				fileNameToChromIDs.get(fileName).add(usedChroms[i]%fileSize);	
			}
		}
		return fileNameToChromIDs;	
	}
	
	//get snp folder
	public static File returnSNPDir(int snp, double maf){
     	String dir_2nd_name = "SNP_" + snp + "_" + Math.floor(maf*1000)/1000;
      	if(dir_2nd_name.charAt(dir_2nd_name.length()-1)=='0')//if the directory name ends with "0.0010", change it to "0.001"
      		dir_2nd_name = dir_2nd_name.substring(0, dir_2nd_name.length()-1);
		File dir_2nd = new File(ParamParser.dataOutputDir, dir_2nd_name);
		return dir_2nd;
	}
	
	public static Integer[] getSNPRealPosition(Integer[] SNPIndex) {
		Arrays.sort(SNPIndex);
		Integer[] snpPos = new Integer[SNPIndex.length];
		//Read the absolute positions of all SNPs from the "out.snp-1.gz" file
		FileHandler readSNPPos = new FileHandler(ParamParser.dataInputDir, ParamParser.posInput);
		readSNPPos.openForRead();
		String nextLine = null;
		int i = 0;
		for(int recordCt = 0;(nextLine = readSNPPos.nextLine()) != null; recordCt++)
			if(Arrays.binarySearch(SNPIndex, recordCt) >= 0)
				snpPos[i++] = Integer.valueOf(nextLine);
		readSNPPos.closeReader();
		return snpPos;
	}
	
	public static Character[] randomSelectAlleleWithoutReplacement() {
		Character[] alleles = {'A', 'T', 'C', 'G'};
		List<Character> alleleList = new ArrayList<Character>(Arrays.asList(alleles));
		Collections.shuffle(alleleList);
		Character[] selectedAlleles = {alleleList.get(0), alleleList.get(1)};
		return  selectedAlleles;
	}
	
	public static void writeOutObject(Object obj, File objectLog){
		try {
			if(objectLog.exists())
				objectLog.delete();
			FileOutputStream localFileOutputStream = new FileOutputStream(objectLog);
			ObjectOutputStream localObjectOutputStream = new ObjectOutputStream(localFileOutputStream);
			localObjectOutputStream.writeObject(obj);
			localObjectOutputStream.flush();
			localObjectOutputStream.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static Object readObject(File objectLog){
		Object obj = null;
		try {
			if(!objectLog.exists())
				System.out.println(objectLog.getName() + "does not exist.");
			FileInputStream fis = new FileInputStream(objectLog);
			ObjectInputStream ois = new ObjectInputStream(fis);
			obj = ois.readObject();
			ois.close();
		} catch (ClassNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return obj;
	}

	public static void generatePlinkInputFiles(String workDir, Map<Integer, byte[]> chromToMarkerSNPs, Integer[] snpPos, List<Subject> subjectList) {
		File workingDir = new File(workDir);
		if(!workingDir.exists()) workingDir.mkdirs();
		FileHandler writeMAPFile = new FileHandler(workDir, ParamParser.outPrefix + "commVar2.map");
		writeMAPFile.openForWrite(true);
		for(int snp : snpPos)			
			writeMAPFile.writeString(ParamParser.founderChrom + "\trs" + snp + "\t0\t" + snp);
		writeMAPFile.closeWriter();			
  		FileHandler writePEDFile = new FileHandler(workDir, ParamParser.outPrefix + "commVar2.ped");
		writePEDFile.openForWrite(true);
		for(Subject s : subjectList){
			byte[] haplotypeA = chromToMarkerSNPs.get(s.getChromA());
			byte[] haplotypeB = chromToMarkerSNPs.get(s.getChromB());
			String genotype = "";
			for(int j = 0; j < snpPos.length; j++)//ArrayList is a sequential list, so its insertion and retrieval order is the same.
				genotype += "\t" + String.valueOf(haplotypeA[j]) + "\t" + String.valueOf(haplotypeB[j]);//Alleles are represented by 1(mutated) and 2 (original)	
			writePEDFile.writeString(s.getID()+"\t"+s.getID()+"\t0\t0\t2\t2" + genotype);	
		}
		writePEDFile.closeWriter();						
	}
	
	public static void outputFamFile(File famFile, List<Subject> caseList, List<Subject> controlList) {
		if(famFile.exists()) return;
		FileHandler writeFam = new FileHandler(famFile);
		writeFam.openForWrite(true);
		int ct = 0;
		int numOfLines = caseList.size() + controlList.size();
		for(Subject caseSub : caseList){
			if(ct < numOfLines) writeFam.writeString(caseSub.getID()+"\t"+caseSub.getID()+"\t0\t0\t2\t2");
			else writeFam.writeLastLine(caseSub.getID()+"\t"+caseSub.getID()+"\t0\t0\t2\t2");
			ct++;
		}
		for(Subject ctrlSub : controlList){
			if(ct < numOfLines) writeFam.writeString(ctrlSub.getID()+"\t"+ctrlSub.getID()+"\t0\t0\t2\t1");
			else writeFam.writeLastLine(ctrlSub.getID()+"\t"+ctrlSub.getID()+"\t0\t0\t2\t1");
			ct++;
		}
		writeFam.closeWriter();
	}
	
	public static void outputMapFile(File mapFile, List<Integer> curTags, int chromID) {
		Collections.sort(curTags);
		FileHandler handleTagSNPs = new FileHandler(mapFile);
		handleTagSNPs.openForWrite(true);
		int ct = 0;
		int numOfLines = curTags.size();
	    for(Integer curTag : curTags){
	    	if(ct < numOfLines) handleTagSNPs.writeString(chromID + "\t" + "rs" + String.valueOf(curTag) + "\t0\t" + String.valueOf(curTag));
	    	else handleTagSNPs.writeLastLine(chromID + "\t" + "rs" + String.valueOf(curTag) + "\t0\t" + String.valueOf(curTag));
	    }
	    handleTagSNPs.closeWriter(); 
	}
}
