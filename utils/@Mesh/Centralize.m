function [K,V2V,E2E] = Centralize(G,scale)
%Centrializes G
%   scale: scale G to a unit 'ScaleArea'

if iscell(G.F)
    error('Not implemented for non-triangular meshes yet');
end
G.V=G.V-repmat(mean(G.V,2),1,size(G.V,2));

if strcmp(scale,'ScaleArea')
    area=G.ComputeSurfaceArea;
    G.V=G.V*sqrt(1/area);
end