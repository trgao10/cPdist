function [TextureCoords1,ImprMap] = TPSDeformation(GM,GN,InputMap,FeatureType,TextureCoords1,TextureCoords2)
%TPSDEFORMATION
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

TPS_FEATURESM = DISCtoPLANE(TextureCoords1(:,TPS_FEATURESM_INDS)','d2p');
TPS_FEATURESN = DISCtoPLANE(TextureCoords2(:,TPS_FEATURESN_INDS)','d2p');
if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
    tP = DISCtoPLANE(TextureCoords1','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    TextureCoords1 = DISCtoPLANE(pt,'p2d')';
    TextureCoords1(:,isnan(compl(TextureCoords1))) = 1;
elseif (length(TPS_FEATURESM)==3) % Affine Transformation
    tP = DISCtoPLANE(TextureCoords1','d2p');
    [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
    pt = [A,b]*[tP';ones(1,size(tP,1))];
    TextureCoords1 = DISCtoPLANE(pt','p2d')';
    TextureCoords1(:,isnan(compl(TextureCoords1))) = 1;
end
ImprMap = knnsearch(TextureCoords2',TextureCoords1');

end

