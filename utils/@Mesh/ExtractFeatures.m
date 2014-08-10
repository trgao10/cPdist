function [Cgauss,Cmean,Cmin,Cmax] = ExtractFeatures(GM,options)
%EXTRACTFEATURES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 2
    options = struct();
end
ConfMaxLocalWidth = getoptions(options,'ConfMaxLocalWidth',8);
GaussMaxLocalWidth = getoptions(options,'GaussMaxLocalWidth',10);
GaussMinLocalWidth = getoptions(options,'GaussMinLocalWidth',6);
ADMaxLocalWidth = getoptions(options,'ADMaxLocalWidth',7);
ExcludeBoundary = getoptions(options,'ExcludeBoundary',1);
Display = getoptions(options,'Display','off');

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% extract features (local maximum of conformal factors)
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
[~,~,Cmin,Cmax,Cmean,Cgauss,~] = GM.ComputeCurvature(options);
DNE = Cmin.^2+Cmax.^2;
DNETruncInds = find(DNE>median(DNE));
[GaussMaxInds,~] = GM.FindLocalMax(Cgauss,GaussMaxLocalWidth,ExcludeBoundary);
[GaussMinInds,~] = GM.FindLocalMax(-Cgauss,GaussMinLocalWidth,ExcludeBoundary);

BB = GM.GetBoundingBox;
diam = sqrt(sum(abs(BB(:,2)-BB(:,1)).^2));
r = 0.3*diam;
AD = AreaDistortionFeature(GM,r);
[ADMaxInds,~] = GM.FindLocalMax(1./AD',ADMaxLocalWidth,1);
ADMaxInds(AD(ADMaxInds)>1) = [];

Ring = GM.ComputeVertexRing;

if isfield(GM.Aux,'Conf')
    [ConfMaxInds,~] = GM.FindLocalMax(GM.Aux.Conf,ConfMaxLocalWidth,ExcludeBoundary);
    ToBeDelInds = [];
    for j=1:length(ConfMaxInds)
        if isempty(find(DNETruncInds == ConfMaxInds(j), 1))
            RingNBD = Ring{ConfMaxInds(j)};
            flag = 0;
            for k=1:length(RingNBD)
                if find(DNETruncInds == RingNBD(k))
                    flag = 1;
                    break;
                end
            end
            if (flag == 0)
                ToBeDelInds = [ToBeDelInds,j];
            end
        end
    end
    ConfMaxInds(ToBeDelInds) = [];
    
    GM.Aux.ConfMaxInds = ConfMaxInds;
end

ToBeDelInds = [];
for j=1:length(GaussMaxInds)
    if isempty(find(DNETruncInds == GaussMaxInds(j), 1))
        ToBeDelInds = [ToBeDelInds,j];
    end
end
GaussMaxInds(ToBeDelInds) = [];

ToBeDelInds = [];
for j=1:length(GaussMinInds)
    if isempty(find(DNETruncInds == GaussMinInds(j), 1))
        ToBeDelInds = [ToBeDelInds,j];
    end
end
GaussMinInds(ToBeDelInds) = [];

ToBeDelInds = [];
for j=1:length(ADMaxInds)
    if isempty(find(DNETruncInds == ADMaxInds(j), 1))
        RingNBD = Ring{ADMaxInds(j)};
        flag = 0;
        for k=1:length(RingNBD)
            if find(DNETruncInds == RingNBD(k))
                flag = 1;
                break;
            end
        end
        if (flag == 0)
            ToBeDelInds = [ToBeDelInds,j];
        end
    end
end
ADMaxInds(ToBeDelInds) = [];

if strcmpi(Display, 'on')
    if isfield(GM.Aux,'name')
        figure('Name',['Features on ' GM.Aux.name]);
    else
        figure;
    end
    GM.draw();hold on;
    set(gcf,'ToolBar','none');
    if isfield(GM.Aux,'Conf')
        scatter3(GM.V(1,ConfMaxInds),GM.V(2,ConfMaxInds),GM.V(3,ConfMaxInds),'g','filled');
    end
    scatter3(GM.V(1,GaussMaxInds),GM.V(2,GaussMaxInds),GM.V(3,GaussMaxInds),'r','filled');
    scatter3(GM.V(1,GaussMinInds),GM.V(2,GaussMinInds),GM.V(3,GaussMinInds),'b','filled');
    scatter3(GM.V(1,ADMaxInds),GM.V(2,ADMaxInds),GM.V(3,ADMaxInds),'y','filled');
    
    set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
    set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
    set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
    set(gca, 'CameraViewAngle', 10.5477);
end

GM.Aux.ADMaxInds = ADMaxInds;
GM.Aux.GaussMaxInds = GaussMaxInds;
GM.Aux.GaussMinInds = GaussMinInds;

end

