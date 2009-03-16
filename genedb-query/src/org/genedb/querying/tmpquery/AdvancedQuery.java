package org.genedb.querying.tmpquery;

import org.genedb.querying.core.QueryClass;
import org.genedb.querying.core.QueryParam;

import org.apache.lucene.index.Term;
import org.apache.lucene.search.BooleanClause;
import org.apache.lucene.search.BooleanQuery;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.TermQuery;
import org.apache.lucene.search.WildcardQuery;
import org.apache.lucene.search.BooleanClause.Occur;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

@QueryClass(
        title="Coding and pseudogenes by protein length",
        shortDesc="Get a list of transcripts ",
        longDesc=""
    )
public class AdvancedQuery extends OrganismLuceneQuery {

    @QueryParam(
            order=1,
            title="Minimum length of protein in bases"
    )
    private String category;

    @QueryParam(
            order=2,
            title="Minimum length of protein in bases"
    )
    private String search = "";

    @Override
    protected String getluceneIndexName() {
        return "org.gmod.schema.mapped.Feature";
    }


    @Override
    protected void getQueryTermsWithoutOrganisms(List<Query> queries) {

        BooleanQuery advQuery = new BooleanQuery();
        if (search.contains("*")) {
            advQuery.add(new WildcardQuery(new Term(category, search.toLowerCase())),Occur.SHOULD);
        } else {
            advQuery.add(new TermQuery(new Term(category, search.toLowerCase())), Occur.SHOULD);
        }

        BooleanQuery booleanQuery = new BooleanQuery();
        booleanQuery.add(new BooleanClause(geneOrPseudogeneQuery, Occur.MUST));
        booleanQuery.add(new BooleanClause(advQuery, Occur.MUST));

        // Organism stuff

        queries.add(booleanQuery);
    }

    @Override
    public Map<String, Object> prepareModelData() {
        Map<String, String> typeMap = new HashMap<String, String>();

        typeMap.put("allNames", "All Names");
        typeMap.put("product", "Products");
        typeMap.put("curatedAnnotation", "Curation/Notes");
        typeMap.put("goTermId", "GO Terms");
        typeMap.put("ecNumber", "EC Number");
        typeMap.put("pfamId", "Pfam");

        Map<String, Object> ret = new HashMap<String, Object>();
        ret.put("typeMap", typeMap);
        return ret;
    }

    // ------ Autogenerated code below here

    public void setSearch(String search) {
        this.search = search;
    }

    public String getSearch() {
        return search;
    }

    public String getCategory() {
        return category;
    }


    public void setCategory(String category) {
        this.category = category;
    }


    @Override
    protected String[] getParamNames() {
        return new String[] {"search"};
    }

}
