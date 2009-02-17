package org.genedb.db.loading;

import org.gmod.schema.feature.Chromosome;

import org.apache.log4j.Logger;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.springframework.context.ApplicationContext;
import org.springframework.context.support.ClassPathXmlApplicationContext;
import org.springframework.transaction.annotation.Transactional;

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;

/**
 * Test the loading of orthologue data in implicit-cluster mode.
 * Implicit cluster mode is engaged when the input file does not
 * contain explicit clusters, but an algorithm is specified (indicating
 * that the orthologues are algorithmically predicted).
 * <p>
 * This is not a unit test, in that it relies on the EMBL loader
 * to load the genes before we load their orthologue data.
 *
 * @author rh11
 *
 */
public class OrthologueLoaderImplicitClusterTest {

    private static final Logger logger = TestLogger.getLogger(OrthologueLoaderImplicitClusterTest.class);

    private static ApplicationContext applicationContext;
    private static OrthologueTester tester;

    private static final String program = "fasta";
    private static final String programVersion = "3.4t26";
    private static final String algorithm = "Reciprocal best match";

    @BeforeClass
    public static void setup() throws IOException, ParsingException {
        applicationContext = new ClassPathXmlApplicationContext(new String[] {"Load.xml", "Test.xml"});

        loadEmblFile("test/data/MRSA252_subset.embl", "Saureus_MRSA252");
        loadEmblFile("test/data/MSSA476_subset.embl", "Saureus_MSSA476");
        loadEmblFile("test/data/EMRSA15_subset.embl", "Saureus_EMRSA15");

        loadOrthologues("test/data/Saureus_subset_genenames.ortho",
            program, programVersion, algorithm, true);
        loadOrthologues("test/data/Saureus_subset_transcriptnames.ortho",
            program, programVersion, algorithm, false);

        tester = applicationContext.getBean("orthologueTester", OrthologueTester.class);
    }

    private static void loadOrthologues(String filename,
            String program, String programVersion, String algorithm,
            boolean geneNames)
        throws IOException, ParsingException {

        OrthologuesLoader loader = applicationContext.getBean("orthologuesLoader", OrthologuesLoader.class);

        File file = new File(filename);
        Reader reader = new FileReader(file);
        try {
            OrthologueFile orthologueFile = new OrthologueFile(file, reader);

            loader.setAnalysisProperties(program, programVersion, algorithm);
            loader.setGeneNames(geneNames);
            loader.load(orthologueFile);
        } finally {
            reader.close();
        }
    }

    @AfterClass
    public static void cleanUp() {
        if (tester == null) {
            // This can happen if there's an error in setup:
            // JUnit still calls us even if setup threw an exception.
            logger.error("Tester is null in cleanUp");
        } else {
            tester.cleanUp();
        }
    }

    private static void loadEmblFile(String filename, String organismCommonName) throws IOException, ParsingException {
        logger.trace(String.format("Loading '%s' into organism '%s'", filename, organismCommonName));

        EmblLoader emblLoader = applicationContext.getBean("emblLoader", EmblLoader.class);
        emblLoader.setOrganismCommonName(organismCommonName);
        emblLoader.setSloppyControlledCuration(true);
        emblLoader.setTopLevelFeatureClass(Chromosome.class);

        File file = new File(filename);
        Reader reader = new FileReader(file);
        try {
            emblLoader.load(new EmblFile(file, reader));
        } finally {
            reader.close();
        }
    }

    private void testOrthologueGroup(Double identity, String... polypeptideUniqueNames) {
        tester.orthologueGroup(program, programVersion, algorithm, identity, polypeptideUniqueNames);
    }

    @Transactional @Test
    public void testGeneNameOrthologueGroups() {
        testOrthologueGroup(100.0, "SAEMRSA1513290:pep", "SAR1478:pep");
        testOrthologueGroup(92.9, "SAEMRSA1516500:pep", "SAR1820:pep");
        testOrthologueGroup(100.0, "SAEMRSA1519750:pep", "SAR2153:pep");
        testOrthologueGroup(95.7, "SAEMRSA1519820:pep", "SAR2160:pep");
        testOrthologueGroup(null, "SAEMRSA1521630:pep", "SAR2349:pep");
        testOrthologueGroup(99.2, "SAEMRSA1502320:pep", "SAR0270:pep");
        testOrthologueGroup(96.6, "SAEMRSA1523410:pep", "SAR2531:pep");
        testOrthologueGroup(100.0, "SAEMRSA1525060:pep", "SAR2680:pep");
        testOrthologueGroup(96.5, "SAEMRSA1503370:pep", "SAR0403:pep");
        testOrthologueGroup(null, "SAEMRSA1504360:pep", "SAR0511:pep");
        testOrthologueGroup(100.0, "SAEMRSA1507570:pep", "SAR0889:pep");
        testOrthologueGroup(100.0, "SAEMRSA1511870:pep", "SAS1280:pep");
        testOrthologueGroup(99.5, "SAEMRSA1517490:pep", "SAS1765:pep");
        testOrthologueGroup(99.6, "SAEMRSA1518070:pep", "SAS1823:pep");
        testOrthologueGroup(null, "SAEMRSA1520150:pep", "SAS2010:pep");
        testOrthologueGroup(98.5, "SAEMRSA1523990:pep", "SAS2388:pep");
        testOrthologueGroup(98.6, "SAEMRSA1524970:pep", "SAS2480:pep");
        testOrthologueGroup(99.2, "SAEMRSA1503330:pep", "SAS0357:pep");
        testOrthologueGroup(98.9, "SAEMRSA1504320:pep", "SAS0463:pep");
        testOrthologueGroup(100.0, "SAEMRSA1504870:pep", "SAS0518:pep");
        testOrthologueGroup(98.6, "SAEMRSA1505550:pep", "SAS0595:pep");
        testOrthologueGroup(98.7, "SAEMRSA1500750:pep", "SAS0083:pep");
        testOrthologueGroup(93.0, "SAEMRSA1508090:pep", "SAS0850:pep");
        testOrthologueGroup(99.4, "SAR1647:pep", "SAS1508:pep");
        testOrthologueGroup(100.0, "SAR1712:pseudogenic_transcript:pep", "SAS1568:pep");
        testOrthologueGroup(99.4, "SAR2389:pep", "SAS2196:pep");
        testOrthologueGroup(100.0, "SAR2601:pep", "SAS2406:pep");
        testOrthologueGroup(97.1, "SAR1812:pep", "SAS1660:pep");
        testOrthologueGroup(99.0, "SAR0736:pep", "SAS0648:pep");
        testOrthologueGroup(99.2, "SAR0156:pep", "SAS0129:pep");
        testOrthologueGroup(97.9, "SAR1639:pep", "SAS1500:pep");
        testOrthologueGroup(100.0, "SAR0015:pep", "SAS0015:pep");
        testOrthologueGroup(100.0, "SAR1663:pep", "SAS1523:pep");
        testOrthologueGroup(100.0, "SAR1939:pep", "SAS1769:pep");
    }

    @Transactional @Test
    public void testPepNameOrthologueGroups() {
        testOrthologueGroup(53.0, "SAEMRSA1511480:pep", "SAR1311:pep");
        testOrthologueGroup(96.9, "SAEMRSA1512940:pep", "SAR1444:pep");
        testOrthologueGroup(99.4, "SAEMRSA1513440:pep", "SAR1493:pep");
        testOrthologueGroup(100.0, "SAEMRSA1514830:pep", "SAR1640:pep");
        testOrthologueGroup(null, "SAEMRSA1514980:pep", "SAR1655:pep");
        testOrthologueGroup(99.6, "SAEMRSA1515150:pep", "SAR1673:pep");
        testOrthologueGroup(85.5, "SAEMRSA1516810:pep", "SAR0692:pep");
        testOrthologueGroup(96.2, "SAEMRSA1516860:pep", "SAR1859:pep");
        testOrthologueGroup(98.5, "SAEMRSA1517760:pep", "SAR1959:pep");
        testOrthologueGroup(null, "SAEMRSA1518420:pep", "SAR2018:pep");
        testOrthologueGroup(99.7, "SAEMRSA1520260:pep", "SAR2206:pep");
        testOrthologueGroup(98.5, "SAEMRSA1521860:pep", "SAR2373:pep");
        testOrthologueGroup(98.5, "SAEMRSA1525490:pep", "SAR2723:pep");
        testOrthologueGroup(100.0, "SAEMRSA1508900:pep", "SAR1032:pep");
        testOrthologueGroup(null, "SAEMRSA1509280:pep", "SAR1072:pep");
        testOrthologueGroup(98.9, "SAEMRSA1510800:pep", "SAS1181:pep");
        testOrthologueGroup(98.9, "SAEMRSA1514330:pep", "SAS1451:pep");
        testOrthologueGroup(99.5, "SAEMRSA1501520:pep", "SAS0162:pep");
        testOrthologueGroup(100.0, "SAEMRSA1519310:pep", "SAS1928:pep");
        testOrthologueGroup(99.7, "SAEMRSA1520710:pep", "SAS2067:pep");
        testOrthologueGroup(95.0, "SAEMRSA1522400:pep", "SAS2234a:pep");
        testOrthologueGroup(99.1, "SAEMRSA1522570:pep", "SAS2250:pep");
        testOrthologueGroup(88.3, "SAEMRSA1523030:pep", "SAS2295:pep");
        testOrthologueGroup(100.0, "SAEMRSA1523490:pep", "SAS2341:pep");
        testOrthologueGroup(99.5, "SAEMRSA1525780:pep", "SAS2557:pep");
        testOrthologueGroup(98.4, "SAEMRSA1504600:pep", "SAS0491:pep");
        testOrthologueGroup(100.0, "SAEMRSA1504710:pep", "SAS0502:pep");
        testOrthologueGroup(98.5, "SAEMRSA1505740:pep", "SAS0613:pep");
        testOrthologueGroup(100.0, "SAR1035:pep", "SAS0997:pep");
        testOrthologueGroup(88.4, "SAR1141:pep", "SAS1101:pep");
        testOrthologueGroup(100.0, "SAR1187:pep", "SAS1145:pep");
        testOrthologueGroup(99.8, "SAR0014:pep", "SAS0014:pep");
        testOrthologueGroup(100.0, "SAR1512:pep", "SAS0939:pep");
        testOrthologueGroup(100.0, "SAR1729:pep", "SAS1585:pep");
        testOrthologueGroup(89.6, "SAR0628:pep", "SAS0587:pep");
        testOrthologueGroup(100.0, "SAR0772:pep", "SAS0684:pep");
        testOrthologueGroup(99.4, "SAR0883:pep", "SAS0791:pep");
        testOrthologueGroup(99.8, "SAR2691:pep", "SAS2498:pep");
        testOrthologueGroup(99.2, "SAR0864:pep", "SAS0773:pep");
        testOrthologueGroup(99.4, "SAR0234:pep", "SAS0217:pep");
        testOrthologueGroup(99.4, "SAR2454:pep", "SAS2256:pep");
        testOrthologueGroup(98.8, "SAR0594:pep", "SAS0547:pep");
    }
}
