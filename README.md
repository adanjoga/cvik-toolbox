**CVIK** is a Cluster Validity Index Tooolbox for automatic determination of the number of clusters (CVIK). It contains more than 70 functions (m-files). This toolbox was developed with MATLAB R2020b.

A first comparison study of the CVIs in **CVIK** using a Differential Evolution approach was reported in:
```
A. José-García and W. Gómez-Flores.
A survey of cluster validity indices for automatic data clustering using differential evolution.
The Genetic and Evolutionary Computation Conference (GECCO '21), Lille, France, 2021.
DOI: 10.1145/3449639.3459341
```

For more information about the used datasets, visit the [mvc-repository](https://mvc-repository.github.io/) website.

---
 
   Cluster validity indices (28)
   -----------------------------
       chindex         - Calinski-Harabasz index (ch).
       cindex          - C index (cind).
       copindex        - COP index (cop).
       csindex         - CS index (cs).
       cvddindex       - Index based on density-involved distance (cvdd).
       cvnnindex       - Index based on nearest neighbors (cvnn).
       dbindex         - Davies-Bouldin index (db).
       db2index        - Enhanced Davies-Bouldin index (db2).
       dbcvindex       - Density-based index (dbcv).
       dunnindex       - Dunn index (dunn).
       gd31index       - Dunn index variant 3,1 (gd31).
       gd33index       - Dunn index variant 3,3 (gd33).
       gd41index       - Dunn index variant 4,1 (gd41).
       gd43index       - Dunn index variant 4,3 (gd43).
       gd51index       - Dunn index variant 5,1 (gd51).
       gd53index       - Dunn index variant 5,3 (gd53).
       lccvindex       - Index based on local cores (lccv).
       pbmindex        - PBM index (pbm).
       sdbwindex       - S_Dbw validity index (sdbw).
       sfindex         - Score Function index (sf).
       silindex        - Silhouette index (sil).
       ssddindex       - Index based on shapes, sizes, densities, and separation distances (ssdd).
       svindex         - SV index (sv).
       symindex        - Symmetry index (sym).
       symdbindex      - Davies-Bouldin index based on symmetry (sdb).
       symdunnindex    - Dunn index based on symmetry (sdi).
       wbindex         - WB index (wb).
       xbindex         - Xie-Beni index (xb).

       cviconfig       - CVI configuration function.
       evalcvi         - CVI evaluation function.


   Proximity (8)
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


   Clustering (4)
   ------------------
       kmedoids        - K-medoids clustering algorithm.
       acde            - An automatic clustering algorithm based on differential evolution.
       tgca            - A two-stage genetic clustering algorithm.
       depso           - An automatic clustering algorithm based on particle swarm optimization.


   Eternal validation (14)
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


## Contact:

```
Adán José-García (adan.jose@cinvestav.mx)
Wilfrido Gómez-Flores (wgomez@cinvestav.mx)
```
