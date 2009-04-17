package org.genedb.web.mvc.model;

import java.util.ArrayList;
import java.util.List;

import junit.framework.Assert;

import org.genedb.db.audit.MockChangeSetImpl;
import org.gmod.schema.feature.Gene;
import org.gmod.schema.feature.Polypeptide;
import org.gmod.schema.feature.Transcript;
import org.junit.Test;
import org.junit.runner.RunWith;
import org.springframework.beans.factory.annotation.Autowired;
import org.springframework.test.context.ContextConfiguration;
import org.springframework.test.context.junit4.SpringJUnit4ClassRunner;


/**
 * To see the final results of the method under test, simply set the cacheSynchroniser.setNoPrintResult(false); 
 * @author larry@sangerinstitute
 *
 */
@RunWith(SpringJUnit4ClassRunner.class)
@ContextConfiguration
public class PeriodicUpdaterTest extends AbstractUpdaterTest{

    @Autowired
    PeriodicUpdater periodicUpdater;
    
    public void testTopLevelFeatureChangeSet(){
        
    }

    
    /**
     * 1.Start by clearing all the caches
     * 2.Initialise the ChangeSet with the polypeptide Ids to be used with various tests
     * @throws Exception
     */
    @Test
    public void testPolypeptideChangeSet()throws Exception{
        Integer newPolyPep = 810;//PFA0010c:pep
        Integer changedPep = 614;//PFA0005w:pep
        
        //Clear all the caches
        CacheSynchroniser cacheSynchroniser = (CacheSynchroniser)periodicUpdater.getIndexUpdaters().get(1);
        cacheSynchroniser.getBmf().getDtoMap().clear();
        cacheSynchroniser.getBmf().getContextMapMap().clear();
        
        //prevent excessive log printing
        cacheSynchroniser.setNoPrintResult(true);
        
        //Get the changeset
        MockChangeSetImpl changeSet = 
            (MockChangeSetImpl)periodicUpdater.getChangeTracker().changes(PeriodicUpdaterTest.class.getName());
        
        //Add new Gene feature to change set
        List<Integer> newFeatureIds = new ArrayList<Integer>();
        changeSet.getNewMap().put(Polypeptide.class, newFeatureIds); 
        newFeatureIds.add(newPolyPep);
        
        //Change Gene feature
        List<Integer> changedFeatureIds = new ArrayList<Integer>();
        changeSet.getChangedMap().put(Polypeptide.class, changedFeatureIds); 
        newFeatureIds.add(changedPep);
        
        
        //execute class under test
        boolean noErrors = periodicUpdater.processChangeSet();
        
        //Assert Transcript DTO with featureID 807 is not null
        //ID 807(PFA0010c:mRNA) is a transcript of Polypeptide 810(PFA0010c:pep)
        TranscriptDTO transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(807);
        Assert.assertNotNull(transcriptDTO);
        Assert.assertEquals("PFA0010c:mRNA", transcriptDTO.getUniqueName());
        
        //Assert Transcript DTO with featureID 611 is not null
        //ID 611(PFA0005w:mRNA) is a transcript of Polypeptide 614(PFA0005w:pep)
        transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(611);
        Assert.assertNotNull(transcriptDTO);
        Assert.assertEquals("PFA0005w:mRNA", transcriptDTO.getUniqueName());
        
        //Assert No severe errors found
        Assert.assertTrue(noErrors);
    }
    
    /**
     * 1.Start by clearing all the caches
     * 2.Initialise the ChangeSet with the transcript Ids to be used with various tests
     * @throws Exception
     */
    @Test
    public void testGeneChangeSet()throws Exception{
        Integer newGeneId = 2;//PFA0170c
        Integer changedGeneId = 610;//PFA0005w
        
        //Clear all the caches
        CacheSynchroniser cacheSynchroniser = (CacheSynchroniser)periodicUpdater.getIndexUpdaters().get(1);
        cacheSynchroniser.getBmf().getDtoMap().clear();
        cacheSynchroniser.getBmf().getContextMapMap().clear();
        
        //prevent excessive log printing
        cacheSynchroniser.setNoPrintResult(true);
        
        //Get the changeset
        MockChangeSetImpl changeSet = 
            (MockChangeSetImpl)periodicUpdater.getChangeTracker().changes(PeriodicUpdaterTest.class.getName());
        
        //Add new Gene feature to change set
        List<Integer> newFeatureIds = new ArrayList<Integer>();
        changeSet.getNewMap().put(Gene.class, newFeatureIds); 
        newFeatureIds.add(newGeneId);
        
        //Change Gene feature in change set
        List<Integer> changedFeatureIds = new ArrayList<Integer>();
        changeSet.getChangedMap().put(Gene.class, changedFeatureIds); 
        newFeatureIds.add(changedGeneId);
        
        
        //execute class under test
        boolean noErrors = periodicUpdater.processChangeSet();
        
        //Assert Transcript DTO with featureID 3 is not null
        //ID 3(PFA0170c:mRNA) is a transcript of Gene ID 2(PFA0170c)
        TranscriptDTO transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(3);
        Assert.assertNotNull(transcriptDTO);
        Assert.assertEquals("PFA0170c:mRNA", transcriptDTO.getUniqueName());
        
        //Assert Transcript DTO with featureID 611 is not null
        //ID 611(PFA0005w:mRNA) is a transcript of Gene ID 610(PFA0005w)
        transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(611);
        Assert.assertNotNull(transcriptDTO);
        Assert.assertEquals("PFA0005w:mRNA", transcriptDTO.getUniqueName());
        
        //Assert No severe errors found
        Assert.assertTrue(noErrors);
        
    }
    
    /**
     * 1.Start by clearing all the caches
     * 2.Initialise the ChangeSet with the transcript Ids to be used with various tests
     * @throws Exception
     */
    @Test
    public void testTranscriptChangeSet()throws Exception{
        Integer newTranscriptId = 7;//PFA0315w:mRNA
        Integer changedTranscriptId = 14;//PFA0380w:mRNA
        Integer deletedTranscriptId = 19;//PFA0440w:mRNA
        
        //Clear all the caches
        CacheSynchroniser cacheSynchroniser = (CacheSynchroniser)periodicUpdater.getIndexUpdaters().get(1);
        cacheSynchroniser.getBmf().getDtoMap().clear();
        cacheSynchroniser.getBmf().getContextMapMap().clear();
        
        //prevent excessive log printing
        cacheSynchroniser.setNoPrintResult(true);
        
        //Get the changeset
        MockChangeSetImpl changeSet = 
            (MockChangeSetImpl)periodicUpdater.getChangeTracker().changes(PeriodicUpdaterTest.class.getName());
        
        //Add new Transcript feature to change set
        List<Integer> newFeatureIds = new ArrayList<Integer>();
        changeSet.getNewMap().put(Transcript.class, newFeatureIds); 
        newFeatureIds.add(newTranscriptId);
        
        //Change  transcript feature in change set
        List<Integer> changedFeatureIds = new ArrayList<Integer>();
        changeSet.getChangedMap().put(Transcript.class, changedFeatureIds); 
        changedFeatureIds.add(changedTranscriptId);
        
        //Delete  transcript feature from change set
        cacheSynchroniser.getBmf().getDtoMap().put(19, new TranscriptDTO());
        List<Integer> deletedFeatureIds = new ArrayList<Integer>();
        changeSet.getDeletedMap().put(Transcript.class, deletedFeatureIds); 
        deletedFeatureIds.add(deletedTranscriptId);
        
        
        //execute class under test
        boolean noErrors = periodicUpdater.processChangeSet();
        
        //Transcript DTO with featureID 7 is not null
        TranscriptDTO transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(newTranscriptId);
        Assert.assertNotNull(transcriptDTO);
        Assert.assertEquals("PFA0315w:mRNA", transcriptDTO.getUniqueName());
        
        //The changed transcript also updates/inserts it's corresponding 
        //TopLevelFeature, in the case the ID is 1
        String contextMap = cacheSynchroniser.getBmf().getContextMapMap().get(1);
        Assert.assertNotNull(contextMap);
        
        //Assert Transcript DTO with featureID 14 is not null
        transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(changedTranscriptId);
        Assert.assertNotNull(transcriptDTO);
        Assert.assertEquals("PFA0380w:mRNA", transcriptDTO.getUniqueName());
        
        //Assert Transcript DTO with featureID 19 IS null
        transcriptDTO = cacheSynchroniser.getBmf().getDtoMap().get(deletedTranscriptId);
        Assert.assertNull(transcriptDTO);
        
        //Assert No severe errors found
        Assert.assertTrue(noErrors);
    }
}
