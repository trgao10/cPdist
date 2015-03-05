function cPdist_landmarkfree_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

options.FeatureType = 'ConfMax';
options.NumDensityPnts = 1000;
options.AngleIncrement = 0.01;
options.NumFeatureMatch = 4;
options.GaussMinMatch = 'on';
options.ProgressBar = 'off';

disp(['Comparing ' GM.Aux.name ' vs ' GN.Aux.name '...']);

rslt = GM.ComputeContinuousProcrustes(GN,options);

cPrslt{str2double(TAXAind1),str2double(TAXAind2)} = rslt;
save(rslt_mat,'cPrslt');

disp(['cPdist(' GM.Aux.name ', ' GN.Aux.name ' = ' num2str(rslt.cPdist) '.']);

end

