package org.genedb.querying.tmpquery;

import org.genedb.querying.core.HqlQuery;
import org.genedb.querying.core.QueryClass;
import org.genedb.querying.core.QueryParam;

import org.apache.log4j.Logger;
import org.hibernate.Query;
import org.springframework.validation.Errors;
import org.springframework.validation.Validator;

@QueryClass(
        title="Transcripts by their type",
        shortDesc="Get a list of transcripts by type",
        longDesc=""
    )
public class ControlledCurationQuery extends OrganismHqlQuery {

    private static final Logger logger = Logger.getLogger(ControlledCurationQuery.class);

    @QueryParam(
            order=1,
            title="Name of CV"
    )
    private BrowseCategory cv;


    @QueryParam(
            order=2,
            title="Name of CV term"
    )
    private String cvTermName;

    @Override
    protected String getHql() {
        return "select f.uniqueName from Feature f, FeatureCvTerm fct where fct.feature=f and fct.cvTerm.name=:cvTermName and fct.cvTerm.cv.name like :cvName order by f.organism";
    }

    @Override
    protected String getOrganismHql() {
        return "and f.organism_id in (:organismList)";
    }
    // ------ Autogenerated code below here

    public BrowseCategory getCv() {
        return cv;
    }


    public void setCv(BrowseCategory cv) {
        this.cv = cv;
    }


    public String getCvTermName() {
        return cvTermName;
    }


    public void setCvTermName(String cvTermName) {
        this.cvTermName = cvTermName;
    }

    @Override
    protected String[] getParamNames() {
        return new String[] {"cvName", "cvTermName"};
    }

    @Override
    protected void populateQueryWithParams(Query query) {
        logger.error(String.format("cvName='%s' cvTermName='%s'", cv.getLookupName(), cvTermName));
        query.setString("cvName", cv.getLookupName());
        query.setString("cvTermName", cvTermName);
    }


            @Override
            @SuppressWarnings("unused")
            public void validate(Object target, Errors errors) {
                return;
            }

            @Override
            @SuppressWarnings("unchecked")
            public boolean supports(Class clazz) {
                return ControlledCurationQuery.class.isAssignableFrom(clazz);
            }


}
