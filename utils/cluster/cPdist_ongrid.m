%%% add: compute landmarkMSE to rslt

function cPdist_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

rslt = GM.ComputeContinuousProcrustes(GN,options);
rslt.lkMSE = ;

cPrslt{TAXAind1,TAXAind2} = rslt;
save(rslt_mat,'cPrslt');

