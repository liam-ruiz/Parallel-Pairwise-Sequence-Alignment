package edu.rice.comp322;

import edu.rice.hj.api.HjDataDrivenFuture;
import edu.rice.hj.api.SuspendableException;

import java.util.ArrayList;

import static edu.rice.hj.Module0.*;
import static edu.rice.hj.Module1.*;

/**
 * A scorer that works in parallel.
 */
public class UsefulParScoring extends AbstractDnaScoring {

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
     * The DDF matrix that models dependencies.
     */
    private final ArrayList<ArrayList<HjDataDrivenFuture<Boolean>>> ddfs = new ArrayList<>();

    /**
     * The # of rows in the DDF matrix.
     */
    private final int ddfRows;
    /**
     * The # of cols in the DDF matrix.
     */
    private final int ddfCols;

    /**
     * The side length of each square chunk.
     */
    private final int chunksize;

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

        // init ddf matrix
        chunksize = 100;
        ddfRows = (int)Math.ceil((double )xLength / chunksize);
        ddfCols = (int)Math.ceil((double)yLength / chunksize);
        for (int i = 0; i < ddfRows; i++) {
            ArrayList<HjDataDrivenFuture<Boolean>> row = new ArrayList<>();
            for (int j = 0; j < ddfCols; j++) {
                row.add(newDataDrivenFuture());
            }
            ddfs.add(row);
        }

    }

    /**
     * Here you should provide an efficient parallel implementation of the Smith-Waterman algorithm that demonstrates
     * real execution time speedup.
     * {@inheritDoc}
     */
    public int scoreSequences(final String x, final String y) throws SuspendableException {

        // TODO: implement this in parallel!
        finish(() -> {
            forseq(0, ddfRows - 1, 0, ddfCols - 1, (i, j) -> {
                // build dependencies list for each chunk
                ArrayList<HjDataDrivenFuture<Boolean>> dependencies = new ArrayList<>();
                // check for diagonal chunk
                if (!(i == 0 || j == 0)) {
                    dependencies.add(ddfs.get(i - 1).get(j - 1));
                }
                // check for above chunk
                if (i != 0) {
                    dependencies.add(ddfs.get(i - 1).get(j));
                }
                // check for left chunk
                if (j != 0) {
                    dependencies.add(ddfs.get(i).get(j - 1));
                }
                asyncAwait(dependencies, () -> {
                    // convert DDF indices to SW indices
                    for (int k = i * chunksize + 1; k <= (i + 1) * chunksize; k++) {
                        for (int l = j * chunksize + 1; l <= (j + 1) * chunksize; l++) {
                            // check that indices are inbounds
                            if (k <= xLength && l <= yLength) {
                                // do dp algorithm and update stored matrix
                                final char XChar = x.charAt(k - 1);
                                final char YChar = y.charAt(l - 1);
                                // find the largest point to jump from, and use it
                                final int diagScore = s[k - 1][l - 1] + getScore(charMap(XChar), charMap(YChar));
                                final int topColScore = s[k - 1][l] + getScore(charMap(XChar), 0);
                                final int leftRowScore = s[k][l - 1] + getScore(0, charMap(YChar));
                                s[k][l] = Math.max(diagScore, Math.max(leftRowScore, topColScore));
                            }
                        }
                    }
                    // put in DDF to mark as done
                    ddfs.get(i).get(j).put(true);
                });
            });
        });
        return s[xLength][yLength];

    }

}

