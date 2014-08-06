function [uniV,uniF] = ComputeMidEdgeUniformization(G,options)
%COMPUTEMIDEDGEUNIFORMIZATION Summary of this function goes here
%   Detailed explanation goes here

if nargin<2
    options = [];
end

SmoothCurvatureFields = getoptions(options,'SmoothCurvatureFields',10);
DensityLocalWidth = getoptions(options,'DensityLocalWidth',5);
ExcludeBoundary = getoptions(options,'ExcludeBoundary',1);

%%% compute mid-edge mesh
[mV,mF,M,E2Vmap] = G.ComputeMidEdge;

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 1) decide where to cut the surface (for flattening/uniformization)
disp('find a face to cut the surface for uniformization');
[v_max_V] = CORR_spread_points_euclidean(G.V',[],200);
GeoField = pdist2(G.V(:,v_max_V)',G.V');

medianGeoField = mean(GeoField,2);
[~, minplc] = min(medianGeoField);
cut_vertex = v_max_V(minplc);
cut_face = find(G.F2V(:,cut_vertex),1);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 2) flatten the mesh conformally 
disp('flatten the mid-edge mesh...')
unmV = CORR_map_mesh_to_plane_nonconforming(G.V',G.F',mF',cut_face,M,E2Vmap,G.nE,0);

unmF = mF';
unmF(cut_face,:) = []; %% it is the same face number as the original mesh

center_ind = cut_vertex;
tind = find(G.F2V(:,center_ind),1);%v_max_V1(kk);
center_ind = mF(1,tind);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 3) map domain to disk (add the infinity point back as sample point)
% transfer the indices of center point to the mid-edge
unmV = CORR_transform_to_disk_new_with_dijkstra(unmV,mF,E2Vmap,center_ind);
unmF = mF';

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% 4) map the original mesh to the disk using the mid-edge structure
disp('flatten the ORIGINAL mesh using the mid-edge flattening...')
uniV = CORR_flatten_mesh_by_COT_using_midedge(G.V',G.F',M,mV,unmF,unmV,cut_face);
uniF = G.F';

G.Aux.UniformizationV = [uniV,zeros(G.nV,1)]';

%%% 5) spread density points on mesh
disp('spread density points on mesh...');
Cgauss = G.ExtractFeatures;
for j=1:SmoothCurvatureFields
    CgaussFace = mean(Cgauss(G.F));
    [~,TriAreas] = G.ComputeSurfaceArea;
    VertAreas = TriAreas'*G.F2V;
    WeightMatrix = repmat(TriAreas,1,G.nV).*G.F2V.*repmat(1./VertAreas,G.nF,1);
    Cgauss = CgaussFace*WeightMatrix;
end
[MaxInds,~] = G.FindLocalMax(Cgauss',DensityLocalWidth,ExcludeBoundary);
[MinInds,~] = G.FindLocalMax(-Cgauss',DensityLocalWidth,ExcludeBoundary);
minds = [MaxInds;MinInds];

spread_pnts = CORR_spread_points_euclidean(G.V',minds,1000);
G.Aux.DensityPnts = [minds; spread_pnts];
G.Aux.SampNumb = length(G.Aux.DensityPnts);

[~,~,Q] = G.PerformFastMarching(G.Aux.DensityPnts);
VorAreaM = zeros(size(G.Aux.DensityPnts));
for i=1:length(G.Aux.DensityPnts)
    VorAreaM(i) = sum(G.Aux.VertArea(Q == G.Aux.DensityPnts(i)));
end
G.Aux.VorArea = VorAreaM;

disp('DONE!');

end

