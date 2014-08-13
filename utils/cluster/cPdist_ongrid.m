function cPdist_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2,LandmarksPath)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 1000;
options.AngleIncrement = 0.01;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'off';
options.ProgressBar = 'off';

disp(['Comparing ' GM.Aux.name ' vs ' GN.Aux.name '...']);

rslt = GM.ComputeContinuousProcrustes(GN,options);
lk2 = GN.V(:,GetLandmarks(GN,LandmarksPath));
lk1 = GN.V(:,rslt.cPmap(GetLandmarks(GM,LandmarksPath)));
rslt.lkMSE = mean(sqrt(sum((lk2-lk1).^2)));

cPrslt{str2num(TAXAind1),str2num(TAXAind2)} = rslt;
save(rslt_mat,'cPrslt');

disp([GM.Aux.name ' vs ' GN.Aux.name ' done.']);

end

