package org.genedb.misc.remapper

import groovy.sql.Sql

// Takes the mapping file produced with ExtractInfo2.groovy
// and remaps the annotation on the old sequence to the new
// sequence

public class MultiChromSimpleLocationRemapper {

    def db
    //def organismId = 19 //remove the hardcoded organism id!
    boolean transcriptNameHack = false

    MultiChromSimpleLocationRemapper(boolean transcriptNameHack) {
        this.transcriptNameHack = transcriptNameHack
        db = Sql.newInstance(
            'jdbc:postgresql://pgsrv3:5432/bigtest',
            'pathdb',
            'LongJ!@n',
            'org.postgresql.Driver')

    }


    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Read the data from the mapping file          */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    private loadNewMappings(String mappingName) {
        def unsorted = []
        File file = new File(mappingName)
        file.eachLine({
            def parts = it.split()
	    //print parts
            def g = [:]
            g['type'] = parts[0]
            g['name'] = parts[1]
            g['contig'] = parts[2]
            g['strand'] = parts[3]
            g['fmin'] = parts[4] as int
            g['fmax'] = parts[5] as int
            unsorted.add(g)
        })
        return unsorted.sort({it.contig})
    }


    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Maps existing features in db to the flocs in mapping file, and prints SQL         */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    private void convert(def features, def featureID, def organismId) {

//      Map contigNameId = lookupContigMap(sequenceNames)

        Map dbids = [:]
        
        List rows = db.rows('''
        select
        f1.uniquename, f1.feature_id, f1.organism_id, cvt.name as type, floc.fmin as fmin, floc.fmax as fmax, floc.strand as strand, floc.featureloc_id as flocid
        from feature f1, featureloc floc, cvterm cvt
        where f1.organism_id=?
        and floc.feature_id=f1.feature_id
        and f1.type_id=cvt.cvterm_id'''
        , [organismId as int])

        // All the featurelocs for this organism in the database
        for (def row in rows) {
            dbids.put(row['uniquename'], [fid:row['feature_id'], type:row['type'], fmin:row['fmin'], fmax:row['fmax'], strand:row['strand'], flocid:row['flocid']])
        }

        def remove = []
        List match = []
        List nomatch = []

        features.each({f ->
          if (!remove.contains(f.type)) {
            if (f.name != "null") {
              String matchedName = nameMatch(dbids, f)
              if (matchedName != null) {
                  int strand = 1
                  if (f['strand'].equals('-')) {
                      strand = -1
                  }
                //println "Match found for '${matchedName}' originally '${f.name}'"
                f['origSrcName'] = dbids.get(matchedName)['name']
                def fid = dbids.get(matchedName)['fid']
                  
                f['mapped'] = true
		//def newContigId = contigNameId.get(f['contig']+"__new")
		//if (newContigId == null) {
		//   println("Unable to get a feature id from map for ${f['contig']}__new")
		   
		//}

		// For now we just use the feature_id to locate the featureloc that needs changing since it's often the case that these features have 1 featureloc but this
		// absolutely needs to be changed!

                def cmd = "update featureloc set fmin=${f['fmin']-1}, fmax=${f['fmax']}, strand=${strand}, srcfeature_id=${featureID} where feature_id=${fid}" //and srcfeature_id in (${currentContigList})"
                println cmd + ";"
                //db.executeUpdate(cmd)
                cmd = "update feature set seqlen=${f['fmax']-f['fmin']+1} where feature_id=${fid}"
                println cmd + ";"
                //db.executeUpdate(cmd)
                match.add(matchedName)
              } else {
                        //println "No match found for '${f.name}'"
                        nomatch.add(f.name)
                    }
                } else {
                    println "f.name is null"
                }
            }
        })

        println "Matched:"
        match.each({
            println it
        })

        println "\n\n\n\n-----------------------------------------------------"
        println "Unmatched:"
        nomatch.each({
            println it
        })
    }


    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */
    /* Method to return the feature_id of a feature */
    /* If fatal is true, exit if the feature is not */
    /* found                                        */
    /* ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ */

    private int lookupFeatureId(String name, boolean fatal) {
        def row = db.firstRow("select feature_id from feature where uniquename=${name}")
        if (row == null) {
            System.err.println "Unable to find '${name}' as an uniquename in db"
            if (fatal) {
	       System.exit(3)
	    }
	    return -1;
        }
        return row[0]
    }


    // This method tries to match the feature names since features get renamed manually and via the EMBL loader.
    // So here we try to 'guess' what the old name would have been
    
    private String nameMatch(def dbids, def f) {
    
        if (dbids.containsKey(f.name)) {
            return f.name
        }
        
        // If the name cannot be found, it does some 'guess work'
        
        // Lets define some regular expression patterns
        
        def trna_transcript = /(\S+):tRNA\.(\d)/
        def ncrna_exon_with_digit =  /(\S+):ncRNA:exon:(\d)/
        def ncrna_exon_without_digit = /(\S+):ncRNA:exon/
	def pseudogenic_transcript = /(\S+):pseudogenic_transcript/
        def pseudogenic_peptide = /(\S+):pseudogenic_transcript:pep/
	def utr = /(\S+):(\d)utr/
	def rrna_exon_with_digit =  /(\S+):rRNA:exon:(\d)/
	def rrna_exon_without_digit =  /(\S+):rRNA:exon/
	def repeat = /\S+:repeat:\d+-\d+/
	def motif = /\S+:motif:\d+/
        def peptide_without_number = /(\S+):pep/


       


        
        // ~~ Non coding exons with a number at the end ~~ (Example: PF01TR002:ncRNA:exon:2 will be turned into PF01TR002:exon:2)
        
        def matcher = ( f.name =~ ncrna_exon_with_digit )
        if(matcher.matches()){
            String tmp
            if(matcher[0][2].equals('1')){
                tmp = matcher[0][1] + ":exon"           
            }else{
                tmp = matcher[0][1] + ":exon:" + matcher[0][2]                        
            }
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
       
        }
        
        // ~~ Non coding exons without a number (TODO: Combine this with method above!!)
        
        matcher = ( f.name =~ ncrna_exon_without_digit )
        if(matcher.matches()){
            String tmp = matcher[0][1] + ":exon"                                   
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
       
        }

	 // ~~ Pseudogenic peptides ~~ (Example: PFA0705c:pseudogenic_transcript:pep gets turned in PFA0705c:pep)
	
	matcher = ( f.name =~ pseudogenic_peptide )
        if(matcher.matches()){
            String tmp = matcher[0][1] + ":pep" ;                             
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                tmp = matcher[0][1] + ":1:pep" ;
		if(dbids.containsKey(tmp)){		
                   return tmp;
                }else {
		   return null;
                }
            }
            
	    
       
        }

        // ~~ Peptides without a number 
        
        matcher = ( f.name =~ peptide_without_number )
        if(matcher.matches()){
            
	  
            String tmp = matcher[0][1] + ":1:pep"                                   
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
       
        }
        
        
        // ~~ Pseudogenic transcripts ~~ (Example: PFA0705c:pseudogenic_transcript gets turned into PFA0705c:mRNA)
	
	matcher = ( f.name =~ pseudogenic_transcript )
        if(matcher.matches()){
            String tmp = matcher[0][1] + ":mRNA"                                   
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
       
        }

   
        
	// ~~ Not-so-case-sensitive UTRs (I know, this is soo hacky! nds 15 Dec 2010)
        
        matcher = (f.name =~ utr)
	if(matcher.matches()){
	    String tmp = matcher[0][1] + ":" + matcher[0][2] + "UTR"
	    if (dbids.containsKey(tmp)) {
                return tmp;
            } 
	}


 	 // ~~ rRNA exons with a number at the end ~~ (Example: PF01TR002:rRNA:exon:2 will be turned into PF01TR002:exon:2)
        
        matcher = ( f.name =~ rrna_exon_with_digit )
        if(matcher.matches()){
            String tmp
            if(matcher[0][2].equals('1')){
                tmp = matcher[0][1] + ":exon"           
            }else{
                tmp = matcher[0][1] + ":exon:" + matcher[0][2]                        
            }
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
       
        }
        
        // ~~ rRNA exons without a number (TODO: Combine this with method above!!)
        
        matcher = ( f.name =~ rrna_exon_without_digit )
        if(matcher.matches()){
            String tmp = matcher[0][1] + ":exon"                                   
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
       
        }


        
        
        if (f.name.endsWith(":tRNA")) {
            String tmp = replaceLast(f.name, ":tRNA", "")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        
 
        
        if (f.name.endsWith(":tRNA.1")) {
            String tmp = replaceLast(f.name, ":tRNA.1", ":tRNA")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":tRNA.1:exon:1")) {
            String tmp = replaceLast(f.name, ":tRNA.1:exon:1", ":exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }

        if (f.name.endsWith(":snoRNA")) {
            String tmp = replaceLast(f.name, ":snoRNA", "")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":snoRNA.1")) {
            String tmp = replaceLast(f.name, ":snoRNA.1", ":snoRNA")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":snoRNA.1:exon:1")) {
            String tmp = replaceLast(f.name, ":snoRNA.1:exon:1", ":exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }

        if (f.name.endsWith(":snRNA.1")) {
            String tmp = replaceLast(f.name, ":snRNA.1", ":snRNA")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":snRNA")) {
            String tmp = replaceLast(f.name, ":snRNA", "")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":snRNA.1:exon:1")) {
            String tmp = replaceLast(f.name, ":snRNA.1:exon:1", ":exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }

        if (f.name.endsWith(":rRNA.1")) {
            String tmp = replaceLast(f.name, ":rRNA.1", ":rRNA")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":rRNA")) {
            String tmp = replaceLast(f.name, ":rRNA", "")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":rRNA.1:exon:1")) {
            String tmp = replaceLast(f.name, ":rRNA.1:exon:1", ":exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        
   

        if (f.name.endsWith(":ncRNA.1")) {
            String tmp = replaceLast(f.name, ":ncRNA.1", ":ncRNA")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":ncRNA")) {
            String tmp = replaceLast(f.name, ":ncRNA.1", "")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(":ncRNA.1:exon:1")) {
            String tmp = replaceLast(f.name, ":ncRNA.1:exon:1", ":exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            } else {
                return null
            }
        }
        if (f.name.endsWith(".1")) {
            String tmp = replaceLast(f.name, ".1", ":mRNA")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1")) {
            String tmp = replaceLast(f.name, ".1", ":pseudogenic_transcript")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:pep")) {
            String tmp = replaceLast(f.name, ".1:pep", ":pep")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:pep")) {
            String tmp = replaceLast(f.name, ".1:pep", ":pseudogenic_transcript:pep")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:1")) {
            String tmp = replaceLast(f.name, ".1:exon:1", ":exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:1")) {
            String tmp = replaceLast(f.name, ".1:exon:1", ":pseudogenic_transcript:exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:2")) {
            String tmp = replaceLast(f.name, ".1:exon:2", ":exon:2")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:2")) {
            String tmp = replaceLast(f.name, ".1:exon:2", ":pseudogenic_transcript:exon:1")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:3")) {
            String tmp = replaceLast(f.name, ".1:exon:3", ":exon:3")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:3")) {
            String tmp = replaceLast(f.name, ".1:exon:3", ":pseudogenic_transcript:exon:3")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:4")) {
            String tmp = replaceLast(f.name, ".1:exon:4", ":exon:4")
            if (dbids.containsKey(tmp)) {
                return tmp;
            }
        }
        if (f.name.endsWith(".1:exon:4")) {
            String tmp = replaceLast(f.name, ".1:exon:4", ":pseudogenic_transcript:exon:4")
            if (dbids.containsKey(tmp)) {
                return tmp
            }
        }
        if (f.name.endsWith(".1:exon:5")) {
            String tmp = replaceLast(f.name, ".1:exon:5", ":exon:5")
            if (dbids.containsKey(tmp)) {
                return tmp
            }
        }
        if (f.name.endsWith(".1:exon:5")) {
            String tmp = replaceLast(f.name, ".1:exon:5", ":pseudogenic_transcript:exon:5")
            if (dbids.containsKey(tmp)) {
                return tmp
            }
        }




        return null
    }

    private String replaceLast(String input, String from, String to) {
        if (!input.endsWith(from)) {
            return input
        }
        int index = input.lastIndexOf(from)
        return input.substring(0, index) + to
    }


    // MAIN METHOD

    /* Takes the following arguments: Mapping file, name of new sequence, organism id and db url */
    
    public static void main(String[] args) {

      if (args.length < 4) {
        println "MultiChromSimpleLocationRemapper mappingFile newSequenceName organism_id dburl"
        System.exit(3)
      }

      MultiChromSimpleLocationRemapper app = new MultiChromSimpleLocationRemapper(true)

      def mapping = [:]
      def currentContigList = "";
      boolean first = true

      // Name of the new sequence (should already be loaded into the database)
      String name = args[1];
      String featureId = app.lookupFeatureId(name, true)
      mapping.put(name, featureId)
      

      String organismId = args[2];
      //for (int i in 1 .. args.length-1) {  
	      // println "Storing mapping '${name}' '${lookup}'"
	       //int oldId = app.lookupContigId(args[i], false) 
	       //if (oldId != -1) {
               // if (!first) {
		//          currentContigList += ", "
	        //    }
	        //   currentContigList += oldId
	        //   first = false
           //} 
      //}

      def features = app.loadNewMappings(args[0])

//      def sequences = {}
//      for (int i in 1 .. args.length-1) {
//          sequences.add(args[i])
//      }
//      app.convert(features, sequences, mapping)
      app.convert(features, featureId, organismId)

    }





}
