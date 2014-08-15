function [rslt] = ImproveMap(GM,GN,DistMatrix,MapMatrix,TaxaCode,options)
%IMPROVEMAP
%   rslt.Gname1:            name of the first mesh
%   rslt.Gname2:            name of the second mesh
%   rslt.ImprDist:            continuous Procrustes distance
%   rslt.ImprMap:             optimal map generating cP distance
%   rslt.invImprMap:          inverse of rslt.ImprMap
%   rslt.TextureCoords1:    texture coordinates for the first mesh
%                           (deformed)
%   rslt.TextureCoords2:    textrue coordinates for the second mesh
%                           (not deformed)
%   rslt.ref:               =0 if Improved map is orientation-preserving
%                           =1 if Improved map is orientation-reversing
%
%   Tingran Gao, trgao10@math.duke.edu
%   last modified: 14 Aug 2014
%

if nargin<6
    options = [];
end

ImprType = getoptions(options,'ImprType','MST');
SmoothMap = getoptions(options,'SmoothMap',1);
Iter = getoptions(options,'Iter','off');

if ~isfield(GM.Aux,'name') && ~isfield(GN.Aux,'name')
    error('Either Mesh missing .Aux.name');
end

rslt.Gname1 = GM.Aux.name;
rslt.Gname2 = GN.Aux.name;

TAXAind = cellfun(@(name) find(strcmpi(TaxaCode,name)),{GM.Aux.name,GN.Aux.name});

switch lower(ImprType)
    case 'mst'
        OptimalPath = FindMSTPath(TAXAind(1),TAXAind(2),DistMatrix);
    case 'spt'
    case 'viterbi'
end

rslt.ImprMap = ComposeMapsAlongPath(OptimalPath,MapMatrix);
[rslt.ImprDist,R] = MapToDist(GM.V,GN.V,rslt.ImprMap,GM.Aux.VertArea);
if det(R)>0
    rslt.ref = 0;
else
    rslt.ref = 1;
end
if SmoothMap==0
    %%% return vertex permutation map
    rslt.TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
    rslt.TextureCoords2(:,isnan(compl(rslt.TextureCoords2))) = ones(2,sum(isnan(compl(rslt.TextureCoords2))));
    if rslt.ref==1
        rslt.TextureCoords2(2,:) = -rslt.TextureCoords2(2,:);
    end
    rslt.TextureCoords1 = rslt.TextureCoords2(:,rslt.ImprMap);
else
    %%% project vertex permutation map to a smooth map
    if strcmpi(Iter,'off')
        [rslt.TextureCoords1,rslt.TextureCoords2,rslt.ImprMap] = ProjMoebius(GM,GN,rslt.ImprMap,rslt.ref,options);
    else
        [rslt.ImprDist,rslt.ImprMap,rslt.TextureCoords1,rslt.TextureCoords2] = IterProjMoebius(GM,GN,rslt.ImprMap,rslt.ref,options);
    end
end

end

