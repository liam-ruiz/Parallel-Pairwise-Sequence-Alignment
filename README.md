There are 3 implementations of a parallelized DP pairwise sequence alignment algorithm.

The first assumes there is no cost of overhead in creating asynchronous tasks (envisions ideal parallelism). This is located in IdealParScoring.java

The second tries to achieve the highest speedup over the sequential algorithm, taking into account actual performance with the overhead of creating asynchronous tasks. This is located in UsefulParScoring.java. 

The third tries to achieve the highest speedup over the sequential algorithm with the added restriction that the algorithm must use less than O(n^2) space. This is located in SparseParScoring.java
