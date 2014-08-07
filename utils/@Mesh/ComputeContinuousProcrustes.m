function [rslt] = ComputeContinuousProcrustes(GM,GN,options)
%COMPUTECONTINUOUSPROCRUSTES Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    options = [];
end

%%% useful shortcuts
compl = @(x) x(1,:)+1i*x(2,:);

%%% feature type for matching
FeatureType = getoptions(options,'FeatureType','ConfMax');
NumDensityPnts = getoptions(options,'NumDensityPnts',100);
AngleIncrement = getoptions(options,'AngleIncrement',0.05);
switch FeatureType
    case 'ADMax'
        FeaturesM = GM.Aux.ADMaxInds;
        FeaturesN = GN.Aux.ADMaxInds;
    case 'GaussMax'
        FeaturesM = GM.Aux.GaussMaxInds;
        FeaturesN = GN.Aux.GaussMaxInds;
    case 'GaussMin'
        FeaturesM = GM.Aux.GaussMinInds;
        FeaturesN = GN.Aux.GaussMinInds;
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
    
    for jj=1:length(FeaturesM)
        progressbar(jj,length(FeaturesM),10);
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
                    % Push features on GM to GN by m
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
                    InterpInds1 = knnsearch(V2',[real(pushInterpCoords1);imag(pushInterpCoords1)]');
                    
                    TPS_DISC_VERTICES_FEATURESM = DISCtoPLANE([real(pushInterpCoords1);imag(pushInterpCoords1)]','d2p');
                    TPS_DISC_VERTICES_FEATURESN = DISCtoPLANE([real(InterpCoords2);imag(InterpCoords2)]','d2p');
                    if length(InterpInds1)>=3
                        if (length(InterpInds1)>3) % TPS (Thin Plate Spline)
                            tP = DISCtoPLANE([real(pushSource);imag(pushSource)]','d2p');
                            [ftps] = TEETH_calc_tps(TPS_DISC_VERTICES_FEATURESM,TPS_DISC_VERTICES_FEATURESN-TPS_DISC_VERTICES_FEATURESM);
                            pt = tP + TEETH_eval_tps(ftps,tP);
                            V1 = DISCtoPLANE(pt,'p2d')';
                        elseif (length(InterpInds1)==3) % affine transformation
                            tP = DISCtoPLANE([real(pushSource);imag(pushSource)]','d2p');
                            [A,b] = PlanarThreePtsDeform(TPS_DISC_VERTICES_FEATURESM,TPS_DISC_VERTICES_FEATURESN);
                            pt = [A,b]*[tP';ones(1,size(tP,1))];
                            V1 = DISCtoPLANE(pt','p2d')';
                        end
                        err = MapToDist(GM.V(:,sourceInds),GN.V(:,targetInds),knnsearch(V2',V1'),VorArea);
                    else
                        err = Inf;
                    end
                end
                %%% Record if best so far
                if ~exist('best_err','var')
                    best_err = err;
                else
                    if (err < best_err)
                        best_err = err;
                        ref12 = ref;
                        best_a = a;
                        best_tet = tet;
                        best_TPS_DISC_VERTICES_FEATURESM = TPS_DISC_VERTICES_FEATURESM;
                        best_TPS_DISC_VERTICES_FEATURESN = TPS_DISC_VERTICES_FEATURESN;
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
if length(best_TPS_DISC_VERTICES_FEATURESM)>=3
    if (length(best_TPS_DISC_VERTICES_FEATURESM)>3) % TPS (Thin Plate Spline)
        tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
        [ftps] = TEETH_calc_tps(best_TPS_DISC_VERTICES_FEATURESM,best_TPS_DISC_VERTICES_FEATURESN-best_TPS_DISC_VERTICES_FEATURESM);
        pt = tP + TEETH_eval_tps(ftps,tP);
        TextureCoords1 = DISCtoPLANE(pt,'p2d')';
    elseif (length(best_TPS_DISC_VERTICES_FEATURESM)==3) % Affine Transformation
        tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
        [A,b] = PlanarThreePtsDeform(best_TPS_DISC_VERTICES_FEATURESM,best_TPS_DISC_VERTICES_FEATURESN);
        pt = [A,b]*[tP';ones(1,size(tP,1))];
        TextureCoords1 = DISCtoPLANE(pt','p2d')';
    end
end

cPmap = knnsearch(TextureCoords2',TextureCoords1');
cPdist = MapToDist(GM.V,GN.V,knnsearch(TextureCoords2',TextureCoords1'),GM.Aux.VertArea);

if ref12==1
    TextureCoords1(2,:) = -TextureCoords1(2,:);
    TextureCoords2(2,:) = -TextureCoords2(2,:);
end

rslt.cPdist = cPdist;
rslt.cPmap = cPmap;
rslt.TextureCoords1 = TextureCoords1;
rslt.TextureCoords2 = TextureCoords2;
rslt.ref12 = ref12;

end

