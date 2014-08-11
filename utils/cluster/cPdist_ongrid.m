function cPdist_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2,LandmarksPath)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

rslt = GM.ComputeContinuousProcrustes(GN,options);
lk2 = G2.V(:,GetLandmarks(G2,LandmarksPath));
lk1 = G2.V(:,rslt.cPmap(GetLandmarks(G1,LandmarksPath)));
rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));

cPrslt{TAXAind1,TAXAind2} = rslt;
save(rslt_mat,'cPrslt');

end