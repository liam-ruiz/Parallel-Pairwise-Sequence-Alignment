package edu.rice.comp322;

import edu.rice.hj.api.SuspendableException;

/**
 * A scorer that works in parallel.
 */
public class UsefulParScoring extends AbstractDnaScoring {

    /**
     * <p>main.</p> Takes the names of two files, and in parallel calculates the sequence aligment scores of the two DNA
     * strands that they represent.
     *
     * @param args The names of two files.
     */
    public static void main(final String[] args) throws Exception {
        final ScoringRunner scoringRunner = new ScoringRunner(UsefulParScoring::new);
        scoringRunner.start("UsefulParScoring", args);
    }

    /**
     * Creates a new UsefulParScoring.
     *
     * @param xLength length of the first sequence
     * @param yLength length of the second sequence
     */
    public UsefulParScoring(final int xLength, final int yLength) {
        if (xLength <= 0 || yLength <= 0) {
            throw new IllegalArgumentException("Lengths (" + xLength + ", " + yLength + ") must be positive!");
        }

        // TODO: implement this!

        throw new UnsupportedOperationException("Parallel allocation not implemented yet!");
    }

    /**
     * Here you should provide an efficient parallel implementation of the Smith-Waterman algorithm that demonstrates
     * real execution time speedup.
     * {@inheritDoc}
     */
    public int scoreSequences(final String x, final String y) throws SuspendableException {

        // TODO: implement this in parallel!

        throw new UnsupportedOperationException("Parallel scoring not implemented yet!");
    }

}

