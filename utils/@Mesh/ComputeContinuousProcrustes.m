function [rslt] = ComputeContinuousProcrustes(GM,GN,options)
%COMPUTECONTINUOUSPROCRUSTES Summary of this function goes here
%   rslt.Gname1:            name of the first mesh
%   rslt.Gname2:            name of the second mesh
%   rslt.cPdist:            continuous Procrustes distance
%   rslt.cPmap:             optimal map generating cP distance
%   rslt.TextureCoords1:    texture coordinates for the first mesh
%                           (deformed)
%   rslt.TextureCoords2:    textrue coordinates for the second mesh
%                           (not deformed)
%   rslt.ref:               =1 if cP map is orientation-preserving
%                           =0 otherwise
%
%   Tingran Gao, trgao10@math.duke.edu
%   last modified: 10 Aug 2014
%


if nargin<3
    options = [];
end

ProgressBar = getoptions(options,'ProgressBar','on');

%%% useful shortcuts
compl = @(x) x(1,:)+1i*x(2,:);

%%% feature type for matching
FeatureType = getoptions(options,'FeatureType','ConfMax');
NumDensityPnts = getoptions(options,'NumDensityPnts',100);
AngleIncrement = getoptions(options,'AngleIncrement',0.05);
NumFeatureMatch = getoptions(options,'NumFeatureMatch',3);
GaussMinMatch = getoptions(options,'GaussMinMatch','on');
switch FeatureType
    case 'ADMax'
        FeaturesM = GM.Aux.ADMaxInds;
        FeaturesN = GN.Aux.ADMaxInds;
    case 'GaussMax'
        FeaturesM = GM.Aux.GaussMaxInds;
        FeaturesN = GN.Aux.GaussMaxInds;
    case 'ConfMax'
        FeaturesM = GM.Aux.ConfMaxInds;
        FeaturesN = GN.Aux.ConfMaxInds;
end

FeaturesMCoords = compl(GM.Aux.UniformizationV(:,FeaturesM));
FeaturesNCoords = compl(GN.Aux.UniformizationV(:,FeaturesN));

%%% check for NaN's in the uniformization of GM
sourceInds = GM.Aux.DensityPnts(1:NumDensityPnts);
source = compl(GM.Aux.UniformizationV(:,sourceInds));
delInds = isnan(source);
source(delInds) = [];
sourceInds(delInds) = [];
VorArea = GM.ComputeVoronoiArea(sourceInds);
%%% check for NaN's in the uniformization of GN
targetInds = GN.Aux.DensityPnts(1:NumDensityPnts);
target = compl(GN.Aux.UniformizationV(:,targetInds));
delInds = isnan(target);
target(delInds) = [];
targetInds(delInds) = [];

for ref=0:1
    if ref==1
        local_target = conj(target);
    else
        local_target = target;
    end
    V2 = [real(local_target);imag(local_target)];
    V2_kdtree = kdtree_build(V2');
    
    for jj=1:length(FeaturesM)
        if strcmpi(ProgressBar,'on')
            progressbar(jj,length(FeaturesM),10);
        end
        for kk=1:length(FeaturesN)
            z_0 = FeaturesMCoords(jj);
            w_0 = FeaturesNCoords(kk);
            if ref==1
                w_0 = conj(w_0);
            end
            
            for tet = 0:AngleIncrement:2*pi %traverse angles
                [a] = CORR_evaluate_disc_moebius_from_tet(tet,z_0,w_0);
                if(a*conj(a) > 0.9999)
                    err = Inf;
                else
                    % push features on GM to GN by m
                    m = [exp(1i*tet) -a*exp(1i*tet); -conj(a) 1];%takes z_0 -> w_0
                    pushFeatureM = CORR_apply_moebius_as_matrix(m,FeaturesMCoords);
                    if ref==0
                        HDist = TEETH_compute_pairwise_hyperbolic_distances(pushFeatureM.',FeaturesNCoords.');
                    elseif ref==1
                        HDist = TEETH_compute_pairwise_hyperbolic_distances(pushFeatureM.',FeaturesNCoords');
                    end
                    [~, tind1] = min(HDist,[],2);
                    [~, tind2] = min(HDist,[],1);
                    tind2 = tind2';
                    InterpInds1 = find(tind2(tind1)==(1:size(HDist,1))');
                    InterpInds2 = tind1(InterpInds1);
                    %%% at the moment, InterpInds1, InterpInds2 are indices
                    %%% on FeaturesM, FeaturesN, respectively
                    InterpCoords1 = FeaturesMCoords(InterpInds1);
                    InterpCoords2 = FeaturesNCoords(InterpInds2);
                    if ref==1
                        InterpCoords2 = conj(InterpCoords2);
                    end
                    
                    %%% now turn InterpInds1, InterpInds2 are into indices
                    %%% on GM, GN, respectively
                    %%% both InterpInds1, InterpInds2 are indices on GN
                    pushSource = CORR_apply_moebius_as_matrix(m,source);
                    pushInterpCoords1 = CORR_apply_moebius_as_matrix(m,InterpCoords1);
                    
                    TPS_DISC_VERTICES_FEATURESM = DISCtoPLANE([real(pushInterpCoords1);imag(pushInterpCoords1)]','d2p');
                    TPS_DISC_VERTICES_FEATURESN = DISCtoPLANE([real(InterpCoords2);imag(InterpCoords2)]','d2p');
                    if length(pushInterpCoords1)>=NumFeatureMatch
                        if (length(pushInterpCoords1)>3) % TPS (Thin Plate Spline)
                            tP = DISCtoPLANE([real(pushSource);imag(pushSource)]','d2p');
                            [ftps] = TEETH_calc_tps(TPS_DISC_VERTICES_FEATURESM,TPS_DISC_VERTICES_FEATURESN-TPS_DISC_VERTICES_FEATURESM);
                            pt = tP + TEETH_eval_tps(ftps,tP);
                            V1 = DISCtoPLANE(pt,'p2d')';
                        elseif (length(pushInterpCoords1)==3) % affine transformation
                            tP = DISCtoPLANE([real(pushSource);imag(pushSource)]','d2p');
                            [A,b] = PlanarThreePtsDeform(TPS_DISC_VERTICES_FEATURESM,TPS_DISC_VERTICES_FEATURESN);
                            pt = [A,b]*[tP';ones(1,size(tP,1))];
                            V1 = DISCtoPLANE(pt','p2d')';
                        end
                        err = MapToDist(GM.V(:,sourceInds),GN.V(:,targetInds),kdtree_nearest_neighbor(V2_kdtree,V1'),VorArea);
                    else
                        err = Inf;
                    end
                end
                %%% Record if best so far
                if ~exist('best_err','var')
                    best_err = err;
                    ref12 = ref;
                    best_a = a;
                    best_tet = tet;
                    TPS_FEATURESM = TPS_DISC_VERTICES_FEATURESM;
                    TPS_FEATURESN = TPS_DISC_VERTICES_FEATURESN;
                else
                    if (err < best_err)
                        best_err = err;
                        ref12 = ref;
                        best_a = a;
                        best_tet = tet;
                        TPS_FEATURESM = TPS_DISC_VERTICES_FEATURESM;
                        TPS_FEATURESN = TPS_DISC_VERTICES_FEATURESN;
                    end
                end
            end
        end
    end
end

m = [exp(1i*best_tet) -best_a*exp(1i*best_tet); -conj(best_a) 1];
pushGM = CORR_apply_moebius_as_matrix(m,compl(GM.Aux.UniformizationV));
pushGM(isnan(pushGM)) = 1+1i;
TextureCoords2 = GN.Aux.UniformizationV(1:2,:);
TextureCoords2(:,isnan(compl(TextureCoords2))) = ones(2,sum(isnan(compl(TextureCoords2))));
if ref12==1
    TextureCoords2(2,:) = -TextureCoords2(2,:);
end

if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
    tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
    [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
    pt = tP + TEETH_eval_tps(ftps,tP);
    TextureCoords1 = DISCtoPLANE(pt,'p2d')';
elseif (length(TPS_FEATURESM)==3) % Affine Transformation
    tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
    [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
    pt = [A,b]*[tP';ones(1,size(tP,1))];
    TextureCoords1 = DISCtoPLANE(pt','p2d')';
end

TextureCoords2_kdtree = kdtree_build(TextureCoords2');
cPmap = kdtree_nearest_neighbor(TextureCoords2_kdtree, TextureCoords1');
% TextureCoords2_kdtree = KDTreeSearcher(TextureCoords2');
% [~,map] = nrsearch(TextureCoords2,TextureCoords1,1,0);
% cPmap = cell2mat(map);
% cPmap = TextureCoords2_kdtree.knnsearch(TextureCoords1');

if strcmpi(GaussMinMatch,'on')
    [~,InterpGaussMinInds2,preInterpGaussMinInds1] = FindMutuallyNearestNeighbors(GM,GN,cPmap,'GaussMin');
    TPS_GaussMinCoords1 = DISCtoPLANE([real(pushGM(preInterpGaussMinInds1));imag(pushGM(preInterpGaussMinInds1))]','d2p');
    TPS_GaussMinCoords2 = DISCtoPLANE(TextureCoords2(:,InterpGaussMinInds2)','d2p');
    TPS_FEATURESM = [TPS_FEATURESM;TPS_GaussMinCoords1];
    TPS_FEATURESN = [TPS_FEATURESN;TPS_GaussMinCoords2];
    if (length(TPS_FEATURESM)>3) % TPS (Thin Plate Spline)
        tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
        [ftps] = TEETH_calc_tps(TPS_FEATURESM,TPS_FEATURESN-TPS_FEATURESM);
        pt = tP + TEETH_eval_tps(ftps,tP);
        TextureCoords1 = DISCtoPLANE(pt,'p2d')';
    elseif (length(TPS_FEATURESM)==3) % Affine Transformation
        tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
        [A,b] = PlanarThreePtsDeform(TPS_FEATURESM,TPS_FEATURESN);
        pt = [A,b]*[tP';ones(1,size(tP,1))];
        TextureCoords1 = DISCtoPLANE(pt','p2d')';
    end
    [~,map] = nrsearch(TextureCoords2,TextureCoords1,1,0);
    cPmap = cell2mat(map);
%     cPmap = TextureCoords2_kdtree.knnsearch(TextureCoords1');
end

cPdist = MapToDist(GM.V,GN.V,cPmap,GM.Aux.VertArea);

if ref12==1
    TextureCoords1(2,:) = -TextureCoords1(2,:);
    TextureCoords2(2,:) = -TextureCoords2(2,:);
end

if isfield(GM.Aux,'name') && isfield(GN.Aux,'name')
    rslt.Gname1 = GM.Aux.name;
    rslt.Gname2 = GN.Aux.name;
end
rslt.cPdist = cPdist;
rslt.cPmap = cPmap;
rslt.TextureCoords1 = TextureCoords1;
rslt.TextureCoords2 = TextureCoords2;
rslt.ref = ref12;

end

function [InterpInds1,InterpInds2,preInterpInds1] = FindMutuallyNearestNeighbors(GM,GN,map,Type)

switch Type
    case 'ADMax'
        GM_MaxInds = GM.Aux.ADMaxInds;
        GN_MaxInds = GN.Aux.ADMaxInds;
    case 'ConfMax'
        GM_MaxInds = GM.Aux.ConfMaxInds;
        GN_MaxInds = GN.Aux.ConfMaxInds;
    case 'GaussMax'
        GM_MaxInds = GM.Aux.GaussMaxInds;
        GN_MaxInds = GN.Aux.GaussMaxInds;
    case 'GaussMin'
        GM_MaxInds = GM.Aux.GaussMinInds;
        GN_MaxInds = GN.Aux.GaussMinInds;
end
pfGM_MaxInds = map(GM_MaxInds);

if ~isempty(GM_MaxInds)&&~isempty(GN_MaxInds)
    [~,~,Q] = GN.PerformFastMarching(pfGM_MaxInds);
    GN2pfGM = Q(GN_MaxInds);
    tind1 = zeros(size(GN_MaxInds));
    for j=1:length(tind1)
        tind1(j) = find(pfGM_MaxInds==GN2pfGM(j));
    end
    
    [~,~,Q] = GN.PerformFastMarching(GN_MaxInds);
    pfGM2GN = Q(pfGM_MaxInds);
    tind2 = zeros(size(pfGM_MaxInds));
    for j=1:length(tind2)
        tind2(j) = find(GN_MaxInds==pfGM2GN(j));
    end
    
    InterpMaxInds2 = find(tind2(tind1)==(1:length(tind1))');
    InterpMaxInds1 = tind1(InterpMaxInds2);
    
    InterpInds1 = pfGM_MaxInds(InterpMaxInds1);
    InterpInds2 = GN_MaxInds(InterpMaxInds2);
    preInterpInds1 = GM_MaxInds(InterpMaxInds1);
end

end

