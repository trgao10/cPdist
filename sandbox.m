G = Mesh('off','../DATA/PNAS/meshes/j14_sas.off');
NaNInds = 1780;

G.DeleteVertex(1780);

figure;G.draw();
hold on;
scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
% scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
G.Write('./j14_sas.off','off',[]);

