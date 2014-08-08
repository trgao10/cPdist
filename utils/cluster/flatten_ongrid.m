function flatten_ongrid(mesh_file, sample_file)
%sample_relax_cP:   sample on one surface
%   Detailed explanation goes here

%==========================================================================
% Preprocessing
%==========================================================================

GM = Mesh('off', mesh_file);
revName = strtok(mesh_file(end:-1:1),'/');
GM.Aux.name = strtok(revName(end:-1:1),'_');
GM.Centralize('ScaleArea');
GM.ComputeMidEdgeUniformization; %%% default options

GM.Nf = GM.ComputeFaceNormals;
GM.Nv = GM.F2V'*GM.Nf';
GM.Nv = GM.Nv'*diag(1./sqrt(sum((GM.Nv').^2,1)));

%%% Compute cotangent Laplacian operator.
GM.Aux.LB = GM.ComputeCotanLaplacian;

%%% Save results to a .mat file.
save(sample_file, 'GM');

end

