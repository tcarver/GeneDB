package org.genedb.web.mvc.model;

import java.io.Serializable;
import java.util.Collections;
import java.util.List;
import java.util.Map;

//import org.apache.log4j.Logger;

import com.thoughtworks.xstream.annotations.XStreamAlias;


public class FeatureDTO implements Serializable {

	//private static final transient Logger logger = Logger.getLogger(TranscriptDTO.class);
	private static final long serialVersionUID = 3878466785198622703L;
	private Map<String,Object> algorithmData;
	
	private boolean anAlternateTranscript;
	private List<String> clusterIds;
	private List<String> comments;
	private List<FeatureCvTermDTO> controlledCurations;
	
	@XStreamAlias("dbxrefs")
	private List<DbXRefDTO> dbXRefDTOs;
	
	private String geneName;
	private List<FeatureCvTermDTO> goBiologicalProcesses;
	private List<FeatureCvTermDTO> goCellularComponents;
	private List<FeatureCvTermDTO> goMolecularFunctions;
	
	
	private long lastModified = Long.MIN_VALUE;
	private int max;
	private int min;
	private List<String> notes;
	private List<String> obsoleteNames;
	private String organismCommonName;
	private String organismHtmlShortName;
	private List<String> orthologueNames;
	
	private List<FeatureCvTermDTO> products;
	private String properName;
	private boolean proteinCoding;
	private boolean pseudo;
	private List<String> publications;
	private short strand;
	
	private Map<String, List<String>> synonymsByTypes;
	private String topLevelFeatureDisplayName;
	private int topLevelFeatureLength;
	private String topLevelFeatureType;
	private String topLevelFeatureUniqueName;
	private String typeDescription;
	private String uniqueName;
	private boolean selenoprotein = false;

	public FeatureDTO() {
		super();
	}

	public boolean isSelenoprotein() {
		return selenoprotein;
	}

	public void setSelenoprotein(boolean selenoprotein) {
		this.selenoprotein = selenoprotein;
	}

//	public boolean isProperGeneName() {
//	    //<c:if test="${!empty dto.geneName && dto.geneName != dto.uniqueName}">
//	
//	    if (getGeneName() == null || getGeneName().length() < 1) {
//	        return false;
//	    }
//	    if (getUniqueName().matches(getGeneName()+":\\w+")) {
//	        return false;
//	    }
//	    if (getUniqueName().matches(getGeneName()+"\\.\\w+")) {
//	        return false;
//	    }
//	    return true;
//	}

	public Map<String, Object> getAlgorithmData() {
	    return algorithmData;
	}

	public List<String> getClusterIds() {
	    return listOrEmptyList(clusterIds);
	}

	public List<String> getComments() {
	    return comments;
	}

	public List<FeatureCvTermDTO> getControlledCurations() {
	    return listOrEmptyList(controlledCurations);
	}

	public List<DbXRefDTO> getDbXRefDTOs() {
	    return listOrEmptyList(dbXRefDTOs);
	}

	public String getGeneName() {
	    return geneName;
	}

	public List<FeatureCvTermDTO> getGoBiologicalProcesses() {
	    return listOrEmptyList(goBiologicalProcesses);
	}

	public List<FeatureCvTermDTO> getGoCellularComponents() {
	    return listOrEmptyList(goCellularComponents);
	}

	public List<FeatureCvTermDTO> getGoMolecularFunctions() {
	    return listOrEmptyList(goMolecularFunctions);
	}

	/**
	 * Get a summary of the protein map for this transcript
	 *
	 * @return the details for the protein map, or null if it shouldn't be shown
	 */
//	public ImageMapSummary getIms() {
//	    return ims;
//	}

	public long getLastModified() {
	    return lastModified;
	}

	public String getLocation() {
	    return String.format("%d - %d %s", min, max, (this.strand < 0) ? " (reverse strand)" : "" );
	}

	public int getMax() {
	    return max;
	}

	public int getMin() {
	    return min;
	}

	public List<String> getNotes() {
	    return notes;
	}

	public List<String> getObsoleteNames() {
	    return listOrEmptyList(obsoleteNames);
	}

	public String getOrganismCommonName() {
	    return organismCommonName;
	}

	public String getOrganismHtmlShortName() {
	    return organismHtmlShortName;
	}

	public List<String> getOrthologueNames() {
	    return listOrEmptyList(orthologueNames);
	}



	public List<FeatureCvTermDTO> getProducts() {
	    return listOrEmptyList(products);
	}

	public String getProperName() {
	    return properName;
	}

	public List<String> getPublications() {
	    return listOrEmptyList(publications);
	}

	public void setSynonymsByTypes(Map<String, List<String>> synonymsByTypes) {
		this.synonymsByTypes = synonymsByTypes;
	}

	public Map<String, List<String>> getSynonymsByTypes() {
	   	return synonymsByTypes;
	}

	public String getTopLevelFeatureDisplayName() {
	    return topLevelFeatureDisplayName;
	}

	public int getTopLevelFeatureLength() {
	    return topLevelFeatureLength;
	}

	public String getTopLevelFeatureType() {
	    return topLevelFeatureType;
	}

	public String getTopLevelFeatureUniqueName() {
	    return topLevelFeatureUniqueName;
	}

	public String getTypeDescription() {
	    return typeDescription;
	}

	public String getUniqueName() {
	    return uniqueName;
	}

	public boolean isAnAlternateTranscript() {
	    return anAlternateTranscript;
	}

	public boolean isProteinCoding() {
	    return proteinCoding;
	}

	public boolean isPseudo() {
	    return pseudo;
	}

	public void setAlgorithmData(Map<String, Object> algorithmData) {
	    this.algorithmData = algorithmData;
	}

	public void setAnAlternateTranscript(boolean anAlternateTranscript) {
	    this.anAlternateTranscript = anAlternateTranscript;
	}

	public void setClusterIds(List<String> clusterIds) {
	    this.clusterIds = clusterIds;
	}

	public void setComments(List<String> comments) {
	    this.comments = comments;
	}

	public void setControlledCurations(List<FeatureCvTermDTO> controlledCurations) {
	    this.controlledCurations = controlledCurations;
	}

	public void setDbXRefDTOs(List<DbXRefDTO> dbXRefDTOs) {
	    this.dbXRefDTOs = dbXRefDTOs;
	}

	

	public void setGeneName(String geneName) {
	    this.geneName = geneName;
	}

	public void setGoBiologicalProcesses(List<FeatureCvTermDTO> goBiologicalProcesses) {
	    this.goBiologicalProcesses = goBiologicalProcesses;
	}

	public void setGoCellularComponents(List<FeatureCvTermDTO> goCellularComponents) {
	    this.goCellularComponents = goCellularComponents;
	}

	public void setGoMolecularFunctions(List<FeatureCvTermDTO> goMolecularFunctions) {
	    this.goMolecularFunctions = goMolecularFunctions;
	}

//	public void setIms(ImageMapSummary ims) {
//	    this.ims = ims;
//	}

	public void setLastModified(long lastModified) {
	    this.lastModified = lastModified;
	}

	public void setMax(int max) {
	    this.max = max;
	}

	public void setMin(int min) {
	    this.min = min;
	}

	public void setNotes(List<String> notes) {
	    this.notes = notes;
	}

	public void setObsoleteNames(List<String> obsoleteNames) {
	    this.obsoleteNames = obsoleteNames;
	}

	public void setOrganismCommonName(String organismCommonName) {
	    this.organismCommonName = organismCommonName;
	}

	public void setOrganismHtmlShortName(String organismHtmlShortName) {
	    this.organismHtmlShortName = organismHtmlShortName;
	}

	public void setOrthologueNames(List<String> orthologueNames) {
	    this.orthologueNames = orthologueNames;
	}

	

	public void setProducts(List<FeatureCvTermDTO> products) {
	    this.products = products;
	}

	public void setProperName(String properName) {
	    this.properName = properName;
	}

	public void setProteinCoding(boolean proteinCoding) {
	    this.proteinCoding = proteinCoding;
	}

	public void setPseudo(boolean pseudo) {
	    this.pseudo = pseudo;
	}

	public void setPublications(List<String> publications) {
	    this.publications = publications;
	}

	public void setStrand(short strand) {
	    this.strand = strand;
	}

	public Short getStrand() {
	    return strand;
	}

//	public void setGeneSynonyms(List<String> geneSynonyms) {
//	    this.geneSynonyms = geneSynonyms;
//	}

	public void setTopLevelFeatureDisplayName(String topLevelFeatureDisplayName) {
	    this.topLevelFeatureDisplayName = topLevelFeatureDisplayName;
	}

	public void setTopLevelFeatureLength(int topLevelFeatureLength) {
	    this.topLevelFeatureLength = topLevelFeatureLength;
	}

	public void setTopLevelFeatureType(String topLevelFeatureType) {
	    this.topLevelFeatureType = topLevelFeatureType;
	}

	public void setTopLevelFeatureUniqueName(String topLevelFeatureUniqueName) {
	    this.topLevelFeatureUniqueName = topLevelFeatureUniqueName;
	}

	public void setTypeDescription(String typeDescription) {
	    this.typeDescription = typeDescription;
	}

	public void setUniqueName(String uniqueName) {
	    this.uniqueName = uniqueName;
	}

	private <T> List<T> listOrEmptyList(List<T> list) {
	    if (list == null || list.size() == 0) {
	        return Collections.emptyList();
	    }
	    return list;
	}

}