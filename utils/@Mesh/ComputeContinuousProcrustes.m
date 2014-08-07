function [rslt] = ComputeContinuousProcrustes(GM,GN)
%COMPUTECONTINUOUSPROCRUSTES Summary of this function goes here
%   Detailed explanation goes here

compl = @(x) x(1,:)+1i*x(2,:);

%%% check for NaN's in the uniformization of GM
source = compl(GM.Aux.UniformizationV);
delInds = isnan(source);
source(delInds) = [];
%%% check for NaN's in the uniformization of GN
target = compl(GN.Aux.UniformizationV);
delInds = isnan(target);
target(delInds) = [];

FeaturesM = GM.Aux.ConfMaxInds;
FeaturesN = GN.Aux.ConfMaxInds;

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
            z_0 = compl(GM.Aux.UniformizationV(:,FeaturesM(jj)));
            w_0 = compl(GN.Aux.UniformizationV(:,FeaturesN(kk)));
            
            for tet = 0:0.05:2*pi %traverse angles: sould be in range 0.05-0.1
                [a] = CORR_evaluate_disc_moebius_from_tet(tet,z_0,w_0);
                if(a*conj(a) > 0.9999)
                    err = Inf;
                else
                    % Push GM to GN by m
                    m = [exp(1i*tet) -a*exp(1i*tet); -conj(a) 1];%takes z_0 -> w_0
                    pushFeatureM = CORR_apply_moebius_as_matrix(m,compl(GM.Aux.UniformizationV(:,FeaturesM)));
                    if ref==0
                        HDist = TEETH_compute_pairwise_hyperbolic_distances(pushFeatureM.',compl(GN.Aux.UniformizationV(:,FeaturesN)).');
                    elseif ref==1
                        HDist = TEETH_compute_pairwise_hyperbolic_distances(pushFeatureM.',compl(GN.Aux.UniformizationV(:,FeaturesN))');
                    end
                    [~, tind1] = min(HDist,[],2);
                    [~, tind2] = min(HDist,[],1);
                    tind2 = tind2';
                    InterpInds1 = find(tind2(tind1)==(1:size(HDist,1))');
                    InterpInds2 = tind1(InterpInds1);
                    InterpInds1 = FeaturesM(InterpInds1);
                    InterpInds2 = FeaturesN(InterpInds2);
                    
                    %%% both InterpInds1, InterpInds2 are indices on GN
                    pushGM = CORR_apply_moebius_as_matrix(m,compl(GM.Aux.UniformizationV));
                    InterpInds1 = knnsearch(V2',[real(pushGM(InterpInds1));imag(pushGM(InterpInds1))]');
                    
                    if length(InterpInds1)>=3
                        if (length(InterpInds1)>3) % TPS (Thin Plate Spline)
                            TPS_DISC_VERTICES = DISCtoPLANE(V2','d2p');
                            tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
                            [ftps] = TEETH_calc_tps(TPS_DISC_VERTICES(InterpInds1,:),TPS_DISC_VERTICES(InterpInds2,:)-TPS_DISC_VERTICES(InterpInds1,:));
                            pt = tP + TEETH_eval_tps(ftps,tP);
                            V1 = DISCtoPLANE(pt,'p2d')';
                        elseif (length(InterpInds1)==3) % affine transformation
                            TPS_DISC_VERTICES = DISCtoPLANE(V2','d2p');
                            tP = DISCtoPLANE([real(pushGM);imag(pushGM)]','d2p');
                            [A,b] = PlanarThreePtsDeform(TPS_DISC_VERTICES(InterpInds1,:),TPS_DISC_VERTICES(InterpInds2,:));
                            pt = [A,b]*[tP';ones(1,size(tP,1))];
                            V1 = DISCtoPLANE(pt','p2d')';
                        end
                        err = MapToDist(GM.V,GN.V,knnsearch(V2',V1'),GM.Aux.VertArea);
                    else
                        err = Inf;
                    end
                end
                % Record if best so far
                if ~exist('best_err','var')
                    best_err = err;
                else
                    if (err < best_err)
                        best_err = err;
                        ref12 = ref;
                        %                     best_a = a;
                        %                     best_tet = tet;
                        TextureCoords1 = V1;
                        TextureCoords2 = V2;
                        best_jj = jj;
                        best_kk = kk;
                    end
                end
            end
        end
    end
end

cPdist = best_err;
cPmap = knnsearch(TextureCoords2',TextureCoords1');

if ref12==1
    TextureCoords1(2,:) = -TextureCoords1(2,:);
    TextureCoords2(2,:) = -TextureCoords2(2,:);
end

rslt.cPdist = cPdist;
rslt.cPmap = cPmap;
rslt.TextureCoords1 = TextureCoords1;
rslt.TextureCoords2 = TextureCoords2;
rslt.ref12 = ref12;

% m = [exp(1i*best_tet) -best_a*exp(1i*best_tet); -conj(best_a) 1];
% ts = compl(GM.Aux.UniformizationV);
% pushGM = CORR_apply_moebius_as_matrix(m,ts);
% pushGM(isnan(pushGM)) = 1+1i;
% V1 = [real(pushGM);imag(pushGM)];
% V2 = GN.Aux.UniformizationV(1:2,:);
% V2(:,isnan(compl(V2))) = ones(2,sum(isnan(compl(V2))));
% if best_ref==1
%     V2(2,:) = -V2(2,:);
% end
% 
% figure;GM.draw();hold on;
% scatter3(GM.V(1,FeaturesM(best_jj)),GM.V(2,FeaturesM(best_jj)),GM.V(3,FeaturesM(best_jj)),50,'g','filled');
% figure;GN.draw();hold on;
% scatter3(GN.V(1,FeaturesN(best_kk)),GN.V(2,FeaturesN(best_kk)),GN.V(3,FeaturesN(best_kk)),50,'g','filled');

end

