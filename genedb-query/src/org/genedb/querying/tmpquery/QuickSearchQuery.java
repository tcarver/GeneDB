package org.genedb.querying.tmpquery;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.apache.log4j.Logger;
import org.apache.lucene.document.Document;
import org.apache.lucene.index.CorruptIndexException;
import org.apache.lucene.index.Term;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.Hit;
import org.apache.lucene.search.Hits;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.WildcardQuery;
import org.apache.lucene.search.BooleanClause.Occur;
import org.genedb.querying.core.QueryException;
import org.genedb.querying.core.QueryParam;
import org.springframework.util.StringUtils;
import org.springframework.validation.Errors;

public class QuickSearchQuery extends OrganismLuceneQuery {

    /**
     * 
     */
    private static final long serialVersionUID = -3007330180211992013L;

    private transient Logger logger  = Logger.getLogger(QuickSearchQuery.class);

    private String searchText;

    @QueryParam(
            order=1,
            title="Search gene products?"
    )
    private boolean product;

    @QueryParam(
            order=1,
            title="Search gene names and synonyms?"
    )
    private boolean allNames;

    @QueryParam(
            order=1,
            title="Include pseudogenes"
    )
    private boolean pseudogenes;

    @Override
    protected String getluceneIndexName() {
        return "org.gmod.schema.mapped.Feature";
    }

    @Override
    protected String[] getParamNames() {
        return new String[] {"searchText", "product", "allNames", "pseudogenes"};
    }

    @Override
    protected void getQueryTermsWithoutOrganisms(List<org.apache.lucene.search.Query> queries) {
        BooleanQuery bq = new BooleanQuery();
        bq.add(new WildcardQuery(new Term("allNames", searchText.toLowerCase())), Occur.SHOULD);
        bq.add(new WildcardQuery(new Term("product", searchText.toLowerCase())), Occur.SHOULD);
        queries.add(bq);
        queries.add(geneOrPseudogeneQuery);
        queries.add(isCurrentQuery);        
    }


    @Override
    protected void getQueryTerms(List<Query> queries) {
        getQueryTermsWithoutOrganisms(queries);
    }

    @Override
    public String getParseableDescription() {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public String getQueryDescription() {
        // TODO Auto-generated method stub
        return null;
    }
    
    @SuppressWarnings("unchecked")
    @Override
    public List getResults() throws QueryException {
       throw new QueryException(
               new Exception("Method Not Implemented"));
    }

    /**
     * Get all the results for the quick query 
     * @return
     * @throws QueryException
     */
    public QuickSearchQueryResults getQuickSearchQueryResults()throws QueryException{
        
        QuickSearchQueryResults quickSearchQueryResults = new QuickSearchQueryResults();        
        List<GeneSummary> geneSummaries = quickSearchQueryResults.getResults();
        TreeMap<String, Integer> taxonGroup = quickSearchQueryResults.getTaxonGroup();
        TreeMap<String, Integer> tempTaxonGroup = new TreeMap<String, Integer>();
        try {
            //taxn name
            List<String> currentTaxonNames = null;
            if (taxons!= null && taxons.length>0){
                currentTaxonNames = taxonNodeManager.getNamesListForTaxons(taxons);
            }
            
            Hits hits = lookupInLucene();

            @SuppressWarnings("unchecked")
            Iterator<Hit> it = hits.iterator();
            while (it.hasNext()) {
                Hit hit = it.next();
                Document document = hit.getDocument();
                 
                //Get the current taxon name from document
                String taxonName = document.get("organism.commonName");
                
                boolean isNoTaxonMatch = currentTaxonNames!= null && !currentTaxonNames.contains(taxonName);
                
                if(taxons==null ){
                    //Categorise the taxons into size of hits
                    populateTaxonGroup(taxonGroup, taxonName);
                    
                    populateGeneSummaries(geneSummaries, document);
                    
                }else if (isNoTaxonMatch){
                    populateTaxonGroup(tempTaxonGroup, taxonName);
                    
                }else{ 
                    populateGeneSummaries(geneSummaries, document);
                    
                    //Categrise the taxons into size of hits
                    populateTaxonGroup(taxonGroup, taxonName);            
                }

            }
            Collections.sort(geneSummaries);
            
            //If no matches are found for current taxon, display all other taxons with a match 
            if (geneSummaries.size() == 0 && taxonGroup.size()==0 && tempTaxonGroup.size()>0){
                taxonGroup.putAll(tempTaxonGroup);
            }
            
            if(currentTaxonNames==null && geneSummaries.size() > 0){
                quickSearchQueryResults.setQuickResultType(QuickResultType.ALL_ORGANISMS_IN_ALL_TAXONS);   
                
            }else if(geneSummaries.size() == 1){
                quickSearchQueryResults.setQuickResultType(QuickResultType.SINGLE_RESULT_IN_CURRENT_TAXON);
                
            }else if (geneSummaries.size() > 1){
                quickSearchQueryResults.setQuickResultType(QuickResultType.MULTIPLE_RESULTS_IN_CURRENT_TAXON);
                
            }else{
                quickSearchQueryResults.setQuickResultType(QuickResultType.NO_EXACT_MATCH_IN_CURRENT_TAXON);
            }
            
        } catch (CorruptIndexException exp) {
            throw new QueryException(exp);
        } catch (IOException exp) {
            throw new QueryException(exp);
        }
        return quickSearchQueryResults;
    }
    
    /**
     * Categrise the taxons into size of hits
     * @param taxonGroup
     * @param taxonName
     */
    private void populateTaxonGroup(TreeMap<String, Integer> taxonGroup, String taxonName){
        Integer currentTaxonHitCount = taxonGroup.get(taxonName);
        if (currentTaxonHitCount==null){
            taxonGroup.put(taxonName, 1);
        }else{                        
            taxonGroup.put(taxonName, ++currentTaxonHitCount);
        }
    }
    
    private void populateGeneSummaries(List<GeneSummary> geneSummaries, Document document){        
        logger.debug(StringUtils.collectionToCommaDelimitedString(document.getFields()));
        GeneSummary gs = convertDocumentToReturnType(document);
        geneSummaries.add(gs);        
    }

    @Override
    public Map<String, Object> prepareModelData() {
        return Collections.emptyMap();
    }

    @Override
    public int getOrder() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public void validate(Object arg0, Errors arg1) {
        // TODO Auto-generated method stub

    }
    
    public class QuickSearchQueryResults{
        private List<GeneSummary> results = new ArrayList<GeneSummary>();
        private TreeMap<String, Integer> taxonGroup = new TreeMap<String, Integer>();
        private QuickResultType quickResultType;
        private String singleResultInTaxonGeneId;
        
        public QuickResultType getQuickResultType() {
            return quickResultType;
        }
        public void setQuickResultType(QuickResultType quickResultType) {
            this.quickResultType = quickResultType;
        }
        
        public List<GeneSummary> getResults() {
            return results;
        }
        public void setResults(List<GeneSummary> results) {
            this.results = results;
        }
        public String getSingleResultInTaxonGeneId() {
            return singleResultInTaxonGeneId;
        }
        public void setSingleResultInTaxonGeneId(String singleResultInTaxonGeneId) {
            this.singleResultInTaxonGeneId = singleResultInTaxonGeneId;
        }
        public TreeMap<String, Integer> getTaxonGroup() {
            return taxonGroup;
        }
        public void setTaxonGroup(TreeMap<String, Integer> taxonGroup) {
            this.taxonGroup = taxonGroup;
        }
    }
    
    public enum QuickResultType{
        ALL_ORGANISMS_IN_ALL_TAXONS,
        NO_EXACT_MATCH_IN_CURRENT_TAXON,
        SINGLE_RESULT_IN_CURRENT_TAXON,
        MULTIPLE_RESULTS_IN_CURRENT_TAXON        
    }

    public String getSearchText() {
        return searchText;
    }

    public void setSearchText(String searchText) {
        this.searchText = searchText;
    }

    public boolean isProduct() {
        return product;
    }

    public void setProduct(boolean product) {
        this.product = product;
    }

    public boolean isAllNames() {
        return allNames;
    }

    public void setAllNames(boolean allNames) {
        this.allNames = allNames;
    }

    public boolean isPseudogenes() {
        return pseudogenes;
    }

    public void setPseudogenes(boolean pseudogenes) {
        this.pseudogenes = pseudogenes;
    }

}