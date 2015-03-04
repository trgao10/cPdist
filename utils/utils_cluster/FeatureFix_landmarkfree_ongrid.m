function FeatureFix_landmarkfree_ongrid(G1,G2,rslt_mat,TAXAind1,TAXAind2,options)

GM = load(G1);
GM = GM.G;
GN = load(G2);
GN = GN.G;

load(rslt_mat);

tic;
disp(['Fixing features for ' GM.Aux.name ' vs ' GN.Aux.name '...']);
rslt = FeatureFix(GM,GN,TAXAind1,TAXAind2,options);
Imprrslt{TAXAind1,TAXAind2} = rslt;
save(rslt_mat,'Imprrslt');
disp(['Feature Fixing for ' GM.Aux.name ' vs ' GN.Aux.name ' done.']);
toc;

end

