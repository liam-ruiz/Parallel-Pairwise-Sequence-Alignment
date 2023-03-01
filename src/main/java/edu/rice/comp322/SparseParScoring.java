package edu.rice.comp322;

import edu.rice.hj.api.SuspendableException;

/**
 * A scorer that allocates memory during computation, so that it may compute the scores for two sequences without
 * requiring O(|X||Y|) memory.
 */
public class SparseParScoring extends AbstractDnaScoring {

    /**
     * <p>main.</p> Takes the names of two files, and in parallel calculates the sequence aligment scores of the two DNA
     * strands that they represent.
     *
     * @param args The names of two files.
     */
    public static void main(final String[] args) throws Exception {
        final ScoringRunner scoringRunner = new ScoringRunner(SparseParScoring::new);
        scoringRunner.start("SparseParScoring", args);
    }

    /**
     * Creates a new SparseParScoring.
     *
     * @param xLength length of the first sequence
     * @param yLength length of the second sequence
     */
    public SparseParScoring(final int xLength, final int yLength) {
        if (xLength <= 0 || yLength <= 0) {
            throw new IllegalArgumentException("Lengths (" + xLength + ", " + yLength + ") must be positive!");
        }

        // TODO: implement this sparsely!

        throw new UnsupportedOperationException("Sparse parallel allocation not implemented yet!");
    }

    /**
     * Here you should implement a parallel version of the SW alignment algorithm that can support datasets where the
     * size of the S matrix exceeds the available memory.
     * {@inheritDoc}
     */
    public int scoreSequences(final String x, final String y) throws SuspendableException {

        // TODO: implement this in parallel!

        throw new UnsupportedOperationException("Parallel scoring not implemented yet!");
    }

}

