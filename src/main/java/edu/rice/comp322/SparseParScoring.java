package edu.rice.comp322;

import edu.rice.hj.api.HjDataDrivenFuture;
import edu.rice.hj.api.SuspendableException;

import java.util.ArrayList;

import static edu.rice.hj.Module0.*;
import static edu.rice.hj.Module1.*;

/**
 * A scorer that allocates memory during computation, so that it may compute the scores for two sequences without
 * requiring O(|X||Y|) memory.
 */
public class SparseParScoring extends AbstractDnaScoring {
    /**
     * The length of the first sequence.
     */
    private final int xLength;
    /**
     * The length of the second sequence.
     */
    private final int yLength;
    /**
     * The array to hold the top row of scores
     */
    private final int[] topArr;
    /**
     * The array to hold the left col of scores
     */
    private final int[] leftArr;
    /**
     *  The array to hold the old diagonal of scores
     */
    private int[] diagOld = null;
    /**
     * The array to hold the most recent diagonal of scores
     */
    private int[] diagRecent = null;

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

        this.xLength = xLength;
        this.yLength = yLength;
        //pre allocate the matrix for alignment, dimension+1 for initializations
        this.leftArr = new int[xLength + 1];
        this.topArr = new int[yLength + 1];


        // init col 0
        for (int ii = 1; ii < xLength + 1; ii++) {
            this.leftArr[ii] = getScore(1, 0) * ii;
        }
        // init row 0
        for (int jj = 1; jj < yLength + 1; jj++) {
            this.topArr[jj] = getScore(0, 1) * jj;
        }




        this.leftArr[0] = 0;
        this.topArr[0] = 0;
        //init diagonal
        this.diagRecent = new int[Math.max(xLength, yLength)];
        this.diagOld = new int[Math.max(xLength, yLength)];

    }

    /**
     * Here you should implement a parallel version of the SW alignment algorithm that can support datasets where the
     * size of the S matrix exceeds the available memory.
     * {@inheritDoc}
     */
    public int scoreSequences(final String x, final String y) throws SuspendableException {

        final int xLength = this.xLength;
        final int yLength = this.yLength;

        HjDataDrivenFuture<Integer> ans = newDataDrivenFuture();

        // for all total sum of indices
        forseq(2, yLength + xLength, (l) -> {
        //for (int l = 2; l <= yLength + xLength; l++) {
            // TODO: ALSO HERE
            final boolean beforeBigDiag = l - 1 <= xLength;
            final int[] temp = new int[Math.max(xLength, yLength)];
            finish(()-> {
                forallChunked(1, l, (ii) -> {
                    // the two characters to be compared
                    final int i = l - ii;
                    final int j = ii;
                    final boolean rowInbounds = i >= 1 && i <= xLength;
                    final boolean colInbounds = j >= 1 && j <= yLength;
                    if (rowInbounds && colInbounds) {
                        final char XChar = x.charAt(i - 1);
                        final char YChar = y.charAt(j - 1);

                        // TODO: MOST LIKELY POINT OF FAILURE
                        final int offset = l > xLength + 1 ? l - xLength - 1 : 0;
                        final int numInDiag = ii - offset - 1;

                        int diagScore;
                        int leftRowScore;
                        int topColScore;

                        // find the largest point to jump from, and use it
                        // get left dependency
                        if (j == 1) {
                            leftRowScore = leftArr[i] + getScore(0, charMap(YChar));
                        } else {
                            final int leftIdx = beforeBigDiag ? numInDiag - 1 : numInDiag;
                            leftRowScore = diagRecent[leftIdx] + getScore(0, charMap(YChar));
                        }
                        // get top dependency
                        if (i == 1) {
                            topColScore = topArr[j] + getScore(charMap(XChar), 0);
                        } else {
                            final int topIdx = beforeBigDiag ? numInDiag : numInDiag + 1;
                            topColScore = diagRecent[topIdx] + getScore(charMap(XChar), 0);
                        }

                        // get diagonal dependency
                        if (i == 1 && j == 1) {
                            diagScore = topArr[0] + getScore(charMap(XChar), charMap(YChar));
                        } else if (i == 1) {
                            diagScore = topArr[j - 1] + getScore(charMap(XChar), charMap(YChar));
                        } else if (j == 1) {
                            diagScore = leftArr[i - 1] + getScore(charMap(XChar), charMap(YChar));
                        } else if (l == xLength + 2) {
                            diagScore = diagOld[numInDiag] + getScore(charMap(XChar), charMap(YChar));
                        } else {
                            final int diagIdx = beforeBigDiag ? numInDiag - 1 : numInDiag + 1;
                            diagScore = diagOld[diagIdx] + getScore(charMap(XChar), charMap(YChar));
                        }



                        final int score = Math.max(diagScore, Math.max(leftRowScore, topColScore));
                        // return case
                        if (i == xLength && j == yLength) {
                            ans.put(score);
                        }
                        temp[numInDiag] = score;
                    }



                });
            });
            diagOld = diagRecent;
            diagRecent = temp;






        });

        return ans.get();
    }

}

