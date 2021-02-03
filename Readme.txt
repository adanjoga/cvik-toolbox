This is the Readme file for the DEMO code related to the submitted manuscript:

"A comparative study of different cluster validity indices using differential evolution"
by Adán José-García and Wilfrido Gómez-Flores



About the CVIs
--------------------------------------------------------------------------
A total of 22 clustering validity indices (CVIs) were implemented. Please refer to the references of each CVI for details. The tested CVIs are:

1. The Xie-Beni index (XB)
2. The Calinski-Harabasz index (CH)
3. The Score Function index (SF)
4. The PBM index (PBM)
5. The CS index (CS)
6. The Dunn index variant 3,1 (gD31)
7. The Dunn index variant 4,1 (gD41)
8. The Dunn index variant 5,1 (gD51)
9. The Dunn index variant 3,3 (gD33)
10. The Dunn index variant 4,3 (gD43)
11. The Dunn index variant 5,3 (gD53)
12. The Enhanced Davies-Bouldin index (mDB)
13. The Davies-Bouldin index (DB)
14. The COP index (COP)
15. The Silhouette index (Sil)
16. The Dunn index (DI)
17. The SV index (SV)
18. The Symmetry index (Sym)
19. The Dunn index based on symmetry (SDI)
20. The Davies-Bouldin index based on symmetry (SDB)
21. The S_Dbw validity index (Sdbw)
22. The C index (CI)

Additionally, four datasets have been provided for testing the CVIs.
1. 'Data_4_3'
2. 'Data_5_2'
3. 'Moon'
4. 'Iris'



About the ACDE algorithm
--------------------------------------------------------------------------
Please refer to the following paper for details about the algorithm:
Authors: Swagatam Das and Ajith Abraham
Paper Title: An improved differential evolution algorithm for automatic clustering
Journal: IEEE Transactions on system, man, and cybernetics — PART A: System and Humans, vol. 38, no. 1, 2008, pp. 218 - 237



How to run the program
---------------------------------------------------------------------------
The "demo_acde.m" file has been provided for running the program on MATLAB 2014b. Edit the file to suit your need. By default, the provided demo attempts to use the Sil index (15) and the dataset "Data_4_3" (1).


About the folders and files
--------------------------------------------------------------------------
demo_acde.m: The main program to run the demo
clustering/acde.m: Implementation of the  ACDE algorithm
cvi/...: This folder contains the routines of the evaluated CVIs
datasets/...: This folder contains the mentioned available datasets
util/...: This folder contains three utility functions in order to compute the euclidean distance (eucdist.m), the adjusted Rand index (pairwiseindex.m), and the function to plot the tested dataset (plotclusters.m)


---------------------------------------------------------------------------
Please feel free to send questions, suggestions, bugs, etc. to ajose@tamps.cinvestav.mx and wgomez@tamps.cinvestav.mx

Adán José-García and Wilfrido Gómez-Flores
15th August 2017
---------------------------------------------------------------------------