package org.genedb.hibernate.search;

import org.gmod.schema.cfg.ChadoAnnotationConfiguration;
import org.gmod.schema.feature.AbstractGene;
import org.gmod.schema.feature.Gap;
import org.gmod.schema.feature.Transcript;
import org.gmod.schema.feature.UTR;
import org.gmod.schema.mapped.Feature;

import org.apache.log4j.Logger;
import org.hibernate.CacheMode;
import org.hibernate.Criteria;
import org.hibernate.FlushMode;
import org.hibernate.ScrollMode;
import org.hibernate.ScrollableResults;
import org.hibernate.SessionFactory;
import org.hibernate.Transaction;
import org.hibernate.criterion.Restrictions;
import org.hibernate.search.FullTextSession;
import org.hibernate.search.Search;
import org.hibernate.search.event.FullTextIndexEventListener;
import org.postgresql.ds.PGSimpleDataSource;

import java.io.Console;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.sql.DataSource;

/**
 * Create Lucene indices.
 * <p>
 * The way it works is as follows. Each indexed feature type is treated in turn; currently
 * the classes <code>AbstractGene</code>, <code>Transcript</code>, <code>UTR</code> and <code>Gap</code>
 * are indexed (in that order). For each type, all features of that type are loaded and indexed
 * in batches of 10.
 * If an exception is thrown while indexing a particular feature, the exception is caught and
 * the whole batch will fail.
 * The members of the failed batch are then put into a queue. When all batches of the relevant type
 * have been processed, the queued members of failed batches are indexed individually. If a feature
 * fails this time, that means it cannot be indexed (due to bad data, or a bug in the code).
 * An error is logged.
 *
 * @author rh11
 */
public class IndexChado {
    private static final Logger logger = Logger.getLogger(IndexChado.class);

    private static void die(String message) {
        System.err.println(message);
        System.err.println();
        System.exit(1);
    }

    /**
     * The number of features to be processed in a single batch. If it's set too
     * high, we run out of heap space.
     */
    private static final int BATCH_SIZE = 10;

    /**
     * Which types of feature to index.
     */
    private static final Collection<Class<? extends Feature>> INDEXED_CLASSES
        = new ArrayList<Class<? extends Feature>>();
    static {
        INDEXED_CLASSES.add(AbstractGene.class);
        INDEXED_CLASSES.add(Transcript.class);
        INDEXED_CLASSES.add(UTR.class);
        INDEXED_CLASSES.add(Gap.class);
        // Add feature types here, if a new type of feature should be indexed.
        // Don't forget to update the class doc comment!
    };

    private String databaseUrl;
    private String databaseUsername;
    private String databasePassword;
    private String indexBaseDirectory;

    /**
     * Create a new instance configured with database connection information taken from
     * system properties. If no database password is supplied, the user is prompted for
     * one on the console.
     *
     * @param batchSize
     */
    public static IndexChado configuredWithSystemProperties() {
        String databaseUrl = System.getProperty("database.url");
        if (databaseUrl == null)
            die("The property database.url must be supplied, "
                    + "e.g. -Ddatabase.url=jdbc:postgres://localhost:10101/malaria_workshop");
        String databaseUsername = System.getProperty("database.user");
        if (databaseUsername == null)
            die("The property database.user must be supplied, " + "e.g. -Ddatabase.user=pathdb");

        String databasePassword = System.getProperty("database.password");
        if (databasePassword == null)
            databasePassword = promptForPassword(databaseUrl, databaseUsername);

        String indexBaseDirectory = System.getProperty("index.base");
        if (indexBaseDirectory == null)
            die("The property index.base must be supplied, "
                    + "e.g. -Dindex.base=/software/pathogen/genedb/indexes");

        return new IndexChado(databaseUrl, databaseUsername, databasePassword, indexBaseDirectory);
    }

    private static String promptForPassword(String databaseUrl, String databaseUsername) {
        Console console = System.console();
        if (console == null)
            die("No password has been supplied, and no console found");

        char[] password = null;
        while (password == null)
            password = console.readPassword("Password for %s@%s: ", databaseUsername, databaseUrl);
        return new String(password);
    }

    /**
     * Create a new instance configured with the specified database connection details.
     *
     * @param databaseUrl
     * @param databaseUsername
     * @param databasePassword
     * @param indexBaseDirectory
     */
    private IndexChado(String databaseUrl, String databaseUsername, String databasePassword,
            String indexBaseDirectory) {
        this.databaseUrl = databaseUrl;
        this.databaseUsername = databaseUsername;
        this.databasePassword = databasePassword;
        this.indexBaseDirectory = indexBaseDirectory;
    }

    private static final Pattern POSTGRES_URL_PATTERN = Pattern.compile("jdbc:postgresql://([^:/]*)(?::(\\d+))?/([^/]+)");
    private DataSource getDataSource() {
        PGSimpleDataSource dataSource = new PGSimpleDataSource();
        Matcher urlMatcher = POSTGRES_URL_PATTERN.matcher(databaseUrl);
        if (!urlMatcher.matches())
            throw new RuntimeException(String.format("Malformed PostgreSQL URL '%s'", databaseUrl));

        dataSource.setServerName(urlMatcher.group(1));
        if (urlMatcher.group(2) != null)
            dataSource.setPortNumber(Integer.parseInt(urlMatcher.group(2)));
        dataSource.setDatabaseName(urlMatcher.group(3));

        dataSource.setUser(databaseUsername);
        dataSource.setPassword(databasePassword);

        return dataSource;
    }

    private Map<Integer,SessionFactory> sessionFactoryByBatchSize = new HashMap<Integer,SessionFactory>();

    /**
     * Get a session factory configured with the database connection information
     * for this instance, and the supplied batch size.
     *
     * @param batchSize
     * @return
     */
    private SessionFactory getSessionFactory(int batchSize) {
        if (sessionFactoryByBatchSize.containsKey(batchSize))
            return sessionFactoryByBatchSize.get(batchSize);

        ChadoAnnotationConfiguration cfg = new ChadoAnnotationConfiguration();
        cfg.setDataSource(getDataSource());
        cfg.configure();

        cfg.setProperty("hibernate.dialect", "org.hibernate.dialect.PostgreSQLDialect");
        cfg.setProperty("hibernate.connection.driver_class", "org.postgresql.Driver");
        cfg.setProperty("hibernate.connection.username", databaseUsername);
        cfg.setProperty("hibernate.connection.password", databasePassword);
        cfg.setProperty("hibernate.connection.url", databaseUrl);

        cfg.setProperty("hibernate.search.default.directory_provider",
            "org.hibernate.search.store.FSDirectoryProvider");
        cfg.setProperty("hibernate.search.worker.batch_size", String.valueOf(batchSize));
        cfg.setProperty("hibernate.search.default.indexBase", indexBaseDirectory);

        FullTextIndexEventListener ft = new FullTextIndexEventListener();
        cfg.setListener("post-insert", ft);
        cfg.setListener("post-update", ft);
        cfg.setListener("post-delete", ft);

        SessionFactory sessionFactory = cfg.buildSessionFactory();
        sessionFactoryByBatchSize.put(batchSize, sessionFactory);
        return sessionFactory;
    }

    /**
     * Create a new session, configured with the supplied batch size.
     *
     * @param batchSize
     * @return
     */
    private FullTextSession newSession(int batchSize) {
        SessionFactory sessionFactory = getSessionFactory(batchSize);
        FullTextSession session = Search.createFullTextSession(sessionFactory.openSession());
        session.setFlushMode(FlushMode.MANUAL);
        session.setCacheMode(CacheMode.IGNORE);
        return session;
    }

    /**
     * Index features of the specified class. First of all indexes the features
     * in batches, and then retries the failures one-by-one.
     *
     * @param featureClass
     * @param numBatches
     */
    public void indexFeatures(Class<? extends Feature> featureClass, int numBatches) {
        FullTextSession session = newSession(BATCH_SIZE);
        Transaction transaction = session.beginTransaction();
        Set<Integer> failed = batchIndexFeatures(featureClass, numBatches, session);
        transaction.commit();
        session.close();

        if (failed.size() > 0)
            reindexFailedFeatures(failed);
    }

    /**
     * Attempt to index features in batches. Returns identifiers of the features
     * that failed to be indexed. (An exception processing a feature will cause
     * the whole batch to fail, so it's worth trying to reindex failed features
     * one-by-one.)
     *
     * @param featureClass the class of features to index
     * @param numBatches the number of batches to process. If zero or negative,
     *                process all
     * @param session
     * @return a set of featureIds of the features that failed to be indexed
     */
    private Set<Integer> batchIndexFeatures(Class<? extends Feature> featureClass,
            int numBatches, FullTextSession session) {

        Set<Integer> failedToLoad = new HashSet<Integer>();
        Criteria criteria = session.createCriteria(featureClass);
        criteria.add(Restrictions.eq("obsolete", false)); // Do not index obsolete features
        if (numBatches > 0)
            criteria.setMaxResults(numBatches * BATCH_SIZE);

        ScrollableResults results = criteria.scroll(ScrollMode.FORWARD_ONLY);

        logger.info(String.format("Indexing %s", featureClass));
        int thisBatchCount = 0;
        Set<Integer> thisBatch = new HashSet<Integer>();

        for (int i = 1; results.next(); i++) {
            Feature feature = (Feature) results.get(0);
            thisBatch.add(feature.getFeatureId());

            boolean failed = false;
            try {
                logger.debug(String.format("Indexing '%s' (%s)", feature.getUniqueName(),
                    feature.getClass()));
                session.index(feature);
            } catch (Exception e) {
                logger.warn("Batch failed", e);
                failed = true;
            }

            if (failed || ++thisBatchCount == BATCH_SIZE) {
                logger.info(String.format("Indexed %d of %s", i, featureClass));
                session.clear();
                thisBatchCount = 0;
                if (failed)
                    failedToLoad.addAll(thisBatch);
                thisBatch = new HashSet<Integer>();
            }
        }
        results.close();
        return failedToLoad;
    }

    /**
     * Attempt to index the provided features individually
     * (i.e. in batches of one). Used to reindex failures
     * from a batch indexing run.
     *
     * @param failed a set of features to reindex
     */
    private void reindexFailedFeatures(Set<Integer> failed) {
        logger.info("Attempting to reindex failed features");
        FullTextSession session = newSession(1);
        Transaction transaction = session.beginTransaction();
        for (int featureId : failed) {
            logger.debug(String.format("Attempting to index feature %d", featureId));
            Feature feature = (Feature) session.load(Feature.class, featureId);
            logger.debug(String.format("Loaded feature '%s'", feature.getUniqueName()));
            try {
                session.index(feature);
                logger.debug("Feature successfully indexed");
            } catch (Exception e) {
                logger.error(String.format("Failed to index feature '%s' on the second attempt",
                    feature.getUniqueName()), e);
            }
            session.clear();
        }
        transaction.commit();
        session.close();
    }


    /**
     * If a command-line argument is passed, it should specify a number of batches.
     * For each feature type, just that many batches are indexed. This is useful for
     * rapid testing: for example if you have changed the structure of the index and
     * want to see whether the resulting indices look plausible. It is not useful for
     * anything else, as far as I know.
     *
     * @param args
     */
    public static void main(String[] args) {
        // The numBatches argument is only useful for quick-and-dirty testing
        int numBatches = -1;
        if (args.length == 1)
            numBatches = Integer.parseInt(args[0]);
        else if (args.length != 0)
            throw new IllegalArgumentException("Unexpected command-line arguments");

        IndexChado indexer = IndexChado.configuredWithSystemProperties();
        for (Class<? extends Feature> featureClass: INDEXED_CLASSES) {
            indexer.indexFeatures(featureClass, numBatches);
        }
    }
}
