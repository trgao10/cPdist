function [TextureCoords1,ImprMap] = LaplacianDeformation(GM,GN,InputMap,FeatureType,TextureCoords1,TextureCoords2)
%LAPLACIANDEFORMATION: Laplacian mesh deformation
%   Current design prefers geodesic mutually nearest neighboring features
%   than Euclidean mutually nearest neighboring features: if a feature is
%   assigned to different target features under these two ad-hoc
%   assignments, then the "Euclidean target" gets "wiped out" in the
%   "unique" section.

[~,TPS_FEATURESN,preTPS_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,FeatureType);
[~,TPS_EUC_FEATURESN,preTPS_EUC_FEATURESM] = FindEuclideanMutuallyNearestNeighbors(GM,GN,InputMap,FeatureType);
TPS_FEATURESM_INDS = [preTPS_FEATURESM;preTPS_EUC_FEATURESM];
TPS_FEATURESN_INDS = [TPS_FEATURESN;TPS_EUC_FEATURESN];
[TPS_FEATURESM_INDS,NoRepeatInds] = unique(TPS_FEATURESM_INDS);
TPS_FEATURESN_INDS = TPS_FEATURESN_INDS(NoRepeatInds);
if ~strcmpi(FeatureType,'GaussMin')
    [~,TPS_GAUSSMIN_FEATURESN,preTPS_GAUSSMIN_FEATURESM] = FindMutuallyNearestNeighbors(GM,GN,InputMap,'GaussMin');
end
TPS_FEATURESM_INDS = [TPS_FEATURESM_INDS;preTPS_GAUSSMIN_FEATURESM];
TPS_FEATURESN_INDS = [TPS_FEATURESN_INDS;TPS_GAUSSMIN_FEATURESN];

KM = Mesh('VF',TextureCoords1,GM.F);
BV = KM.FindBoundaries;
anchors = [TPS_FEATURESM_INDS;BV];
L = KM.ComputeCotanLaplacian;
L(anchors,:) = 0;
for j=1:length(anchors)
    L(anchors(j),anchors(j)) = 1;
end
rhs = zeros(KM.nV,2);
rhs(BV,:) = TextureCoords1(:,BV)';
rhs(TPS_FEATURESM_INDS,:) = TextureCoords2(:,TPS_FEATURESN_INDS)';

TextureCoords1 = (L\rhs)';

ImprMap = knnsearch(TextureCoords2',TextureCoords1');

end

