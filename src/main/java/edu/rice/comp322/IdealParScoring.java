package edu.rice.comp322;

import edu.rice.hj.api.SuspendableException;

import static edu.rice.hj.Module1.*;

/**
 * A scorer that works in parallel.
 */
public class IdealParScoring extends AbstractDnaScoring {
    /**
     * The length of the first sequence.
     */
    private final int xLength;
    /**
     * The length of the second sequence.
     */
    private final int yLength;

    /**
     * The Smith-Waterman matrix.
     */
    private final int[][] s;

    /**
     * <p>main.</p> Takes the names of two files, and in parallel calculates the sequence aligment scores of the two DNA
     * strands that they represent.
     *
     * @param args The names of two files.
     */
    public static void main(final String[] args) throws Exception {
        final ScoringRunner scoringRunner = new ScoringRunner(IdealParScoring::new);
        scoringRunner.start("IdealParScoring", args);

    }

    /**
     * Creates a new IdealParScoring.
     *
     * @param xLength length of the first sequence
     * @param yLength length of the second sequence
     */
    public IdealParScoring(final int xLength, final int yLength) throws SuspendableException {
        if (xLength <= 0 || yLength <= 0) {
            throw new IllegalArgumentException("Lengths (" + xLength + ", " + yLength + ") must be positive!");
        }

        this.xLength = xLength;
        this.yLength = yLength;
        //pre allocate the matrix for alignment, dimension+1 for initializations
        this.s = new int[xLength + 1][yLength + 1];


        // init col 0
        for (int ii = 1; ii < xLength + 1; ii++) {
            s[ii][0] = getScore(1, 0) * ii;
        }
        // init row 0
        for (int jj = 1; jj < yLength + 1; jj++) {
            s[0][jj] = getScore(0, 1) * jj;
        }



        //init diagonal
        this.s[0][0] = 0;


    }

    /**
     * This method should be filled in with a parallel implementation of the Smith-Waterman alignment algorithm that
     * maximizes ideal parallelism.
     * {@inheritDoc}
     */
    public int scoreSequences(final String x, final String y) throws SuspendableException {

        // TODO: implement this in parallel!
        final int xLength = this.xLength;
        final int yLength = this.yLength;
        final int[][] matrix = this.s;


        // for all total sum of indices
        forseq(2, yLength + xLength, (l) -> {
            forall(1, l, (ii) -> {
                // the two characters to be compared
                final int i = l - ii;
                final int j = ii;
                final boolean rowInbounds = i >= 1 && i <= xLength;
                final boolean colInbounds = j >= 1 && j <= yLength;
                if (rowInbounds && colInbounds) {

                    final char XChar = x.charAt(i - 1);
                    final char YChar = y.charAt(j - 1);

                    // find the largest point to jump from, and use it
                    final int diagScore = matrix[i - 1][j - 1] + getScore(charMap(XChar), charMap(YChar));
                    final int topColScore = matrix[i - 1][j] + getScore(charMap(XChar), 0);
                    final int leftRowScore = matrix[i][j - 1] + getScore(0, charMap(YChar));
                    doWork(1);
                    matrix[i][j] = Math.max(diagScore, Math.max(leftRowScore, topColScore));
                }

            });

        });
        return matrix[xLength][yLength];

    }

}

