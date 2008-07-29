package org.gmod.schema.feature;

import org.gmod.schema.cfg.FeatureType;
import org.gmod.schema.mapped.Feature;
import org.gmod.schema.mapped.Organism;

import org.biojava.bio.symbol.Location;

import java.sql.Timestamp;

import javax.persistence.Entity;

@Entity
@FeatureType(cv="sequence", term="centromere")
public class Centromere extends Region {

    Centromere() {
        // empty
    }

    public Centromere(Organism organism, String systematicId, boolean analysis,
            boolean obsolete, Timestamp dateAccessioned) {
        super(organism, systematicId, analysis, obsolete, dateAccessioned);
    }

    public static Centromere make(Feature parent, Location exonLocation,
            String exonSystematicId, Organism organism, Timestamp now) {

        Centromere centromere = new Centromere(organism, exonSystematicId, false, false, now);
        return centromere;
    }
}
