function [V1,V2,proj_map12] = ProjMoebius(GM,GN,map12,ref12,options)
%PROJMOEBIUS: Find the Moebius transform from s1 to s2 that is closest to
%             map12. In some sense this is "projecting a map to the
%             subspace of Moebius transforms
%   GM:       triangular mesh with M vertices
%   GN:       triangular mesh with N vertices
%   map12:    Mx1 vector, max(map12)<=N, so that GN.V(:,map12) becomes a
%             vector of size Mx1, which can be used as the texture
%             coordinates for GM
%   ref12:    0 (no reflection) or 1 (reflection)
%
%   Tingran Gao, Duke University
%   trgao10@math.duke.edu
%

compl = @(x) x(1,:)+1i*x(2,:);

if nargin<5
    options = [];
end
GaussMinInds = getoptions(options,'GaussMinInds','off');

%%% check for NaN's in the uniformization of GM
ts = GM.Aux.UniformizationV(1,:)+1i*GM.Aux.UniformizationV(2,:);
delInds = isnan(ts);
ts(delInds) = [];
%%% check for NaN's in the uniformization of GM
V2 = GN.Aux.UniformizationV(1:2,:);
V2(:,isnan(compl(V2))) = ones(2,sum(isnan(compl(V2))));
ref_GM = V2(1:2,map12);
ref_GM(:,delInds) = [];

if ref12==1
    V2(2,:) = -V2(2,:);
    ref_GM(2,:) = -ref_GM(2,:);
end

FeaturesM = GM.Aux.ConfMaxInds;
FeaturesN = GN.Aux.ConfMaxInds;

if options.debug==1
    pfFeaturesM = map12(FeaturesM);
    [D,~,~] = GN.PerformFastMarching(FeaturesN);
    options.method = 'continuous';
    Paths = compute_geodesic_mesh(D, GN.V, GN.F, pfFeaturesM, options);
    
    colmap = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
    colmap = [colmap;colmap*0.6;colmap*0.4;colmap*0.2];
    figure;
    h(1) = subplot(1,2,1);
    GM.draw();
    hold on;
    for j=1:length(FeaturesM)
        scatter3(GM.V(1,FeaturesM(j)),GM.V(2,FeaturesM(j)),GM.V(3,FeaturesM(j)),50,colmap(j,:),'filled');
    end
    h(2) = subplot(1,2,2);
    GN.draw();
    hold on;
    for j=1:length(pfFeaturesM)
        scatter3(GN.V(1,pfFeaturesM(j)),GN.V(2,pfFeaturesM(j)),GN.V(3,pfFeaturesM(j)),50,colmap(j,:),'filled');
    end
    scatter3(GN.V(1,FeaturesN),GN.V(2,FeaturesN),GN.V(3,FeaturesN),30,'w','filled');
    for j=1:length(Paths)
        geoPath = Paths{j};
        plot3(geoPath(1,:),geoPath(2,:),geoPath(3,:),'Color','r','LineWidth',5);
    end
    
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
    set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
    set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
    set(gca, 'CameraViewAngle', 10.5477);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% starts Moebius projection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
best_err = Inf;
best_a = 0;
best_tet = 0;

for jj=1:length(FeaturesM)
    progressbar(jj,length(FeaturesM),10);
    for kk=1:length(FeaturesN)
        z_0 = GM.Aux.UniformizationV(1,FeaturesM(jj))+1i*GM.Aux.UniformizationV(2,FeaturesM(jj));
        w_0 = V2(1,FeaturesN(kk))+1i*V2(2,FeaturesN(kk));
        
        for tet=0:0.01:2*pi %traverse angles: sould be in range 0.05-0.1
            [a] = CORR_evaluate_disc_mobius_from_tet(tet,z_0,w_0);
            if(a*conj(a) > 0.9999)
                err = Inf;
            else
                % Push GM to GN by m
                m = [exp(1i*tet) -a*exp(1i*tet); -conj(a) 1];%takes z_0 -> w_0
                push_GM = CORR_apply_mobius_as_matrix(m,ts);
                push_GM = [real(push_GM);imag(push_GM)];
                err = sum(sum((push_GM-ref_GM).^2));
            end
            % Record if best so far
            if (err < best_err)
                best_err = err;
                best_a = a;
                best_tet = tet;
%                 best_jj = jj;
%                 best_kk = kk;
            end
        end
    end
end

% figure;GM.draw();hold on;
% scatter3(GM.V(1,FeaturesM(best_jj)),GM.V(2,FeaturesM(best_jj)),GM.V(3,FeaturesM(best_jj)),30,'g','filled');
% figure;GN.draw();hold on;
% scatter3(GN.V(1,FeaturesN(best_kk)),GN.V(2,FeaturesN(best_kk)),GN.V(3,FeaturesN(best_kk)),30,'g','filled');

m = [exp(1i*best_tet) -best_a*exp(1i*best_tet); -conj(best_a) 1];
ts = GM.Aux.UniformizationV(1,:)+1i*GM.Aux.UniformizationV(2,:);
push_GM = CORR_apply_mobius_as_matrix(m,ts);
push_GM(isnan(push_GM)) = 1+1i;
V1 = [real(push_GM);imag(push_GM)];

proj_map12 = knnsearch(V2',V1');

if options.debug==1
    pfFeaturesM = proj_map12(FeaturesM);
    options.method = 'continuous';
    Paths = compute_geodesic_mesh(D, GN.V, GN.F, pfFeaturesM, options);
    
    colmap = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
    colmap = [colmap;colmap*0.6;colmap*0.4;colmap*0.2];
    figure;
    h(1) = subplot(1,2,1);
    GM.draw();
    hold on;
    for j=1:length(FeaturesM)
        scatter3(GM.V(1,FeaturesM(j)),GM.V(2,FeaturesM(j)),GM.V(3,FeaturesM(j)),50,colmap(j,:),'filled');
    end
    h(2) = subplot(1,2,2);
    GN.draw();
    hold on;
    for j=1:length(pfFeaturesM)
        scatter3(GN.V(1,pfFeaturesM(j)),GN.V(2,pfFeaturesM(j)),GN.V(3,pfFeaturesM(j)),50,colmap(j,:),'filled');
    end
    scatter3(GN.V(1,FeaturesN),GN.V(2,FeaturesN),GN.V(3,FeaturesN),30,'w','filled');
    for j=1:length(Paths)
        geoPath = Paths{j};
        plot3(geoPath(1,:),geoPath(2,:),geoPath(3,:),'Color','r','LineWidth',5);
    end
    
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
    set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
    set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
    set(gca, 'CameraViewAngle', 10.5477);
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% find mutually nearest maximal Area Distortion points
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[InterpADMaxInds1,InterpADMaxInds2] = FindMutuallyNearestNeighbors(GM,GN,proj_map12,'ADMax');

if options.debug==1
    GM_ADMaxInds = GM.Aux.ADMaxInds;
    pfGM_ADMaxInds = proj_map12(GM_ADMaxInds);
    GN_ADMaxInds = GN.Aux.ADMaxInds;
    [D,~,~] = GN.PerformFastMarching(GN_ADMaxInds);
    options.method = 'continuous';
    Paths = compute_geodesic_mesh(D, GN.V, GN.F, pfGM_ADMaxInds, options);
    
    colmap = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
    colmap = [colmap;colmap*0.6;colmap*0.4;colmap*0.2];
    figure;
    h(1) = subplot(1,2,1);
    GM.draw();
    hold on;
    for j=1:length(GM_ADMaxInds)
        scatter3(GM.V(1,GM_ADMaxInds(j)),GM.V(2,GM_ADMaxInds(j)),GM.V(3,GM_ADMaxInds(j)),50,colmap(j,:),'filled');
    end
    h(2) = subplot(1,2,2);
    GN.draw();
    hold on;
    for j=1:length(pfGM_ADMaxInds)
        scatter3(GN.V(1,pfGM_ADMaxInds(j)),GN.V(2,pfGM_ADMaxInds(j)),GN.V(3,pfGM_ADMaxInds(j)),50,colmap(j,:),'filled');
    end
    scatter3(GN.V(1,GN_ADMaxInds),GN.V(2,GN_ADMaxInds),GN.V(3,GN_ADMaxInds),30,'w','filled');
    for j=1:length(Paths)
        geoPath = Paths{j};
        plot3(geoPath(1,:),geoPath(2,:),geoPath(3,:),'Color','r','LineWidth',5);
    end
    
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
    set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
    set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
    set(gca, 'CameraViewAngle', 10.5477);
    
    figure;
    GN.draw();
    hold on;
    for j=1:length(InterpADMaxInds1)
        scatter3(GN.V(1,InterpADMaxInds1(j)),GN.V(2,InterpADMaxInds1(j)),GN.V(3,InterpADMaxInds1(j)),50,colmap(j,:),'filled');
        scatter3(GN.V(1,InterpADMaxInds2(j)),GN.V(2,InterpADMaxInds2(j)),GN.V(3,InterpADMaxInds2(j)),50,colmap(j,:),'filled');
    end
    
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% find mutually nearest maximal Conformal Factor points
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% find mutually closest minimal Gaussian Curvature points
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if strcmpi(GaussMinInds,'on')
    disp('Matching GaussMinInds.');
    [InterpGaussMinInds1,InterpGaussMinInds2] = FindMutuallyNearestNeighbors(GM,GN,proj_map12,'GaussMin');
    
    %     pfGM_MaxInds = InterpGaussMinInds1;
    %     GN_MaxInds = InterpGaussMinInds2;
    %
    %     [D,~,~] = GN.PerformFastMarching(GN_MaxInds);
    %     options.method = 'continuous';
    %     Paths = compute_geodesic_mesh(D, GN.V, GN.F, pfGM_MaxInds, options);
    %
    %     colmap = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
    %     colmap = [colmap;colmap*0.6;colmap*0.4;colmap*0.2];
    %     h(1) = subplot(1,2,1);
    %     GM.draw();
    %     hold on;
    %     for j=1:length(pfGM_MaxInds)
    %         GM_MaxInds = find(proj_map12==pfGM_MaxInds(j));
    %         scatter3(GM.V(1,GM_MaxInds),GM.V(2,GM_MaxInds),GM.V(3,GM_MaxInds),50,colmap(j,:),'filled');
    %     end
    %     h(2) = subplot(1,2,2);
    %     GN.draw();
    %     hold on;
    %     for j=1:length(pfGM_MaxInds)
    %         scatter3(GN.V(1,pfGM_MaxInds(j)),GN.V(2,pfGM_MaxInds(j)),GN.V(3,pfGM_MaxInds(j)),50,colmap(j,:),'filled');
    %     end
    %     scatter3(GN.V(1,GN_MaxInds),GN.V(2,GN_MaxInds),GN.V(3,GN_MaxInds),30,'w','filled');
    %     for j=1:length(Paths)
    %         geoPath = Paths{j};
    %         plot3(geoPath(1,:),geoPath(2,:),geoPath(3,:),'Color','r','LineWidth',5);
    %     end
    %
    %     Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    %     setappdata(gcf, 'StoreTheLink', Link);
    %
    %     set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
    %     set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
    %     set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
    %     set(gca, 'CameraViewAngle', 10.5477);
end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% TPS: both InterpInds1, InterpInds2 are indices on GN
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if strcmpi(GaussMinInds,'off')
    InterpInds1 = InterpADMaxInds1;
    InterpInds2 = InterpADMaxInds2;
else
    InterpInds1 = [InterpADMaxInds1;InterpGaussMinInds1];
    InterpInds2 = [InterpADMaxInds2;InterpGaussMinInds2];
end

[~,NonRepInds,~] = unique(InterpInds1);
InterpInds1 = InterpInds1(NonRepInds);
InterpInds2 = InterpInds2(NonRepInds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% colmap = [1,0,0;0,1,0;0,0,1;1,1,0;1,0,1;0,1,1];
% colmap = [colmap;colmap*0.6;colmap*0.4;colmap*0.2];
% figure;
% h(1) = subplot(1,2,1);
% GM.draw();
% hold on;
% for j=1:length(InterpInds1)
%     preimg = find(proj_map12==InterpInds1(j));
%     scatter3(GM.V(1,preimg),GM.V(2,preimg),GM.V(3,preimg),50,colmap(j,:),'filled');
% end
% h(2) = subplot(1,2,2);
% GN.draw();
% hold on;
% for j=1:length(InterpInds1)
%     scatter3(GN.V(1,InterpInds1(j)),GN.V(2,InterpInds1(j)),GN.V(3,InterpInds1(j)),50,colmap(j,:),'filled');
%     scatter3(GN.V(1,InterpInds2(j)),GN.V(2,InterpInds2(j)),GN.V(3,InterpInds2(j)),30,colmap(j,:),'filled');
% end
% 
% Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
% setappdata(gcf, 'StoreTheLink', Link);
% 
% set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
% set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
% set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
% set(gca, 'CameraViewAngle', 10.5477);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (length(InterpInds1)>3) % TPS (Thin Plate Spline)
    TPS_DISC_VERTICES = DISCtoPLANE(V2','d2p');
    tP = DISCtoPLANE(V1','d2p');
    [ftps] = TEETH_calc_tps(TPS_DISC_VERTICES(InterpInds1,:),TPS_DISC_VERTICES(InterpInds2,:)-TPS_DISC_VERTICES(InterpInds1,:));
    pt = tP + TEETH_eval_tps(ftps,tP);
    V1 = DISCtoPLANE(pt,'p2d')';
elseif (length(InterpInds1)==3) % affine transformation
    TPS_DISC_VERTICES = DISCtoPLANE(V2','d2p');
    tP = (DISCtoPLANE(V1','d2p'))';
    [A,b] = PlanarThreePtsDeform(TPS_DISC_VERTICES(InterpInds1,:),TPS_DISC_VERTICES(InterpInds2,:));
    pt = [A,b]*[tP;ones(1,size(V1,2))];
    V1 = DISCtoPLANE(pt','p2d')';
else % (length(InterpInds1)<3)
    error('Not Enough Faithful Features to Initiate Deformation!');
end

proj_map12 = knnsearch(V2',V1');

end

function [InterpInds1,InterpInds2] = FindMutuallyNearestNeighbors(GM,GN,map,Type)

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
end

end
