# WORKFLOW DESCRIPTION #

#### CONTINUOUS PROCRUSTES MINIMAL SPANNING TREE (CP-MST) DISTANCE ##

1. Run "cluster_flatten.m" to generate samples (saved under results/).

2. Run "cluster_cPdist.m" to generate pairwise comparision results (saved under rslts/; cluster scripts and debugging information stored under cluster/).

3. Run "cPProcessRslts.m" to separate results from rslts/. The results are stored in results/; typically cPdistMatrix.mat, cPmapsMartix.mat, lmkMSEMatrix.m, TextureCoords1/, TextureCoords2/ are generated and saved.

4. Move rslts/, TextureCoords1/, TextureCoords2/ to other locations (taking too much space).

5. Run "cluster_Imprdist.m" to try different types of improvements.

6. Run "ImprProcessRslts.m" to separate results from rslts/. The results are stored in results/; typically cP[ImprType]distMatrix.mat, cP[ImprType]MapsMartix.mat, cP[ImprType]lmkMSEMatrix.m, TextureCoords1/, TextureCoords2/ are generated and saved.

7. As in Step 4, move rslts/, TextureCoords1/, TextureCoords2/ to other locations for future reference.

---------------------------------
Tingran Gao, Duke University
Aug 23, 2014

