package org.genedb.jogra.services;

import org.genedb.db.taxon.TaxonNode;
import org.genedb.jogra.domain.Term;

import java.sql.SQLException;
import java.util.List;

public interface TermService {
    
    /* Accepts multiple organisms and returns a list of corresponding terms */
    List<Term> getTerms(List<TaxonNode> taxonNode, String termType);
   
    /* Takes the new term, list of old terms and the corrected text (if any) */
    RationaliserResult rationaliseTerm(List<Term> oldTerms, String correctedText, boolean allOrganisms, List<TaxonNode> selectedTaxons) throws SQLException;
       
    /* Method to retrieve systematic IDs (restricted to selected scope) for a given term */
    List<String> getSystematicIDs(Term term, List<TaxonNode> selectedTaxons);
    
    /* Method to retrieve evidence codes for a given term */
    List<String> getEvidenceCodes(Term term);
 
}