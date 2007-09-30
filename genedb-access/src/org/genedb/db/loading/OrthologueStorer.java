package org.genedb.db.loading;

import org.gmod.schema.cv.CvTerm;
import org.gmod.schema.sequence.Feature;
import org.gmod.schema.sequence.FeatureRelationship;

import org.springframework.stereotype.Repository;
import org.springframework.transaction.annotation.Propagation;
import org.springframework.transaction.annotation.Transactional;

import java.io.File;
import java.util.List;
import java.util.Map;


public interface OrthologueStorer {

	void afterPropertiesSet();

	void process(File[] files);

	void writeToDb();

	@Transactional(propagation=Propagation.REQUIRES_NEW)
	public void storePairs(GenePair pair, CvTerm relationship);

	@Transactional(propagation=Propagation.REQUIRES_NEW)
	public void storeCluster(Map.Entry<String, List<String>> entry, CvTerm relationship);
	
	@Transactional(propagation = Propagation.REQUIRES_NEW)
    public void finishClusterHandling();
	
}
