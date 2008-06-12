package org.genedb.web.tags.db;

import static javax.servlet.jsp.PageContext.APPLICATION_SCOPE;
import static org.genedb.web.mvc.controller.TaxonManagerListener.TAXON_NODE_MANAGER;

import org.genedb.db.taxon.TaxonNode;
import org.genedb.db.taxon.TaxonNodeManager;

import java.io.IOException;

import javax.servlet.jsp.JspException;
import javax.servlet.jsp.JspWriter;
import javax.servlet.jsp.tagext.SimpleTagSupport;

public abstract class AbstractHomepageTag extends SimpleTagSupport {

	private String top = "Root"; // FIXME
    
	
    @Override
    public void doTag() throws JspException, IOException {
        TaxonNodeManager tnm = (TaxonNodeManager) 
            getJspContext().getAttribute(TAXON_NODE_MANAGER, APPLICATION_SCOPE);
        
        TaxonNode topNode = tnm.getTaxonNodeForLabel(top);
        
        if (topNode == null) {
            throw new JspException("Homepage Tag: Can't identify taxonNode for '"+top+"'");
        }
        
        JspWriter out = getJspContext().getOut();
        display(topNode, out, 0);
    }
    
    abstract protected void display(TaxonNode node, JspWriter out, int indent) throws IOException;

	public void setTop(String top) {
		this.top = top;
	}

}
