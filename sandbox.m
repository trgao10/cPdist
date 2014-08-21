G = Mesh('off','../DATA/PNAS/meshes/V13_sas.off');
NaNInds = 3969;

G.DeleteVertex(3969);

figure;G.draw();
hold on;
scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
G.Write('./V13_sas.off','off',[]);
