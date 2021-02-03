The Cluster Validity Index Tooolbox (CVIK) for automatic determination of clusters contains more than 70 functions (m-files).
This toolbox was developed with MATLAB R2020b.

Developed by
   + Author 1
   + Author 2

Please, cite the following paper where this toolbox was introduced:

   [PENDING]


**IMPORTANT**: First run the "RUN_ME_FIRST.m" file to add this toolbox to search path. Also, look at the demo file to start using the CVIK toolbox.

---
 
   CVI (Cluster validity indices - 24)
   -----------------------------
       chindex         - Calinski-Harabasz index (ch).
       cindex          - C index (cind).
       copindex        - COP index (cop).
       csindex         - CS index (cs).
       dbindex         - Davies-Bouldin index (db).
       db2index        - Enhanced Davies-Bouldin index (db2).
       dunnindex       - Dunn index (dunn).
       gindex          - Gamma index (gamma).
       gd31index       - Dunn index variant 3,1 (gd31).
       gd33index       - Dunn index variant 3,3 (gd33).
       gd41index       - Dunn index variant 4,1 (gd41).
       gd43index       - Dunn index variant 4,3 (gd43).
       gd51index       - Dunn index variant 5,1 (gd51).
       gd53index       - Dunn index variant 5,3 (gd53).
       osindex         - OS index (os).
       pbmindex        - PBM index (pbm).
       sdbwindex       - S_Dbw validity index (sdbw).
       sfindex         - Score Function index (sf).
       silindex        - Silhouette index (sil).
       svindex         - SV index (sv).
       symindex        - Symmetry index (sym).
       symdbindex      - Davies-Bouldin index based on symmetry (sdb).
       symdunnindex    - Dunn index based on symmetry (sdi).
       xbindex         - Xie-Beni index (xb).

       cviconfig       - CVI configuration function.
       evalcvi         - CVI evaluation function.


   PROXIMITY
   -------------------
       eucdist         - Euclidean distance (euc).
       neucdist        - Normalized Euclidean distance (neuc).
       cosdist         - Cosine similarity (cos).
       pcorr           - Pearson's correlation coefficient (pcorr).
       scorr           - Spearman's correlation coefficient (scorr).
       lapdist         - Laplacian distance (lap).
       symdist         - Symilarity-based distance (sym).
       medist          - Maxium Edge distance (med).

       proxconfig      - Proximity configuration function


   CLUSTERING
   ------------------
       kmedoids        - K-medoids clustering algorithm.
       acde            - An automatic clustering algorithm based on differential evolution.
       tgca            - A two-stage genetic clustering algorithm.
       depso           - An automatic clustering algorithm based on particle swarm optimization.


   EXTERNAL VALIDATION
   -------------------
       inftheoryindex  - External validity indices based on information theory:
                           - Mutual information (mi).
                           - Mutual information (mi).
                           - Variation of information (vi).
                           - Normalized mutual information (nvi).

       pairwiseindex   - External validity indices based on pairwise similarity:
                           - Rand index (ri).
                           - Adjusted rand index (ari).
                           - Wallace coefficient A->B (wab).
                           - Wallace coefficient B->A (wba).
                           - Jaccard index (jrd).
                           - Fowlkes-Mallows index (fm).
                           - Larsen index A->B (lab).
                           - Larsen index B->A (lba).
                           - Meila-Heckerman index (mh).
                           - Mirkin coefficient (mc).

