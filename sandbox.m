G = Mesh('off','../DATA/PNAS/meshes/j17_sas.off');
NaNInds = [3118,3121,3122];

G.DeleteVertex(3121);

figure;G.draw();
hold on;
scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
G.Write('./j17_sas.off','off',[]);