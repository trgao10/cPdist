%%% add: compute landmarkMSE


function [rslt] = cPdist_ongrid(sample_file_1,sample_file_2,rslt_mat_filename)

GM = load(sample_file_1);
GM = GM.G;
GN = load(sample_file_2);
GN = GN.G;

load(rslt_mat_filename);

rslt = GM.ComputeContinuousProcrustes(GN,options);

%output results in the dist matrix and save it back
cPrslt{TAXAind1,TAXAind2} = rslt;
save(rslt_mat_filename,'cPrslt');

