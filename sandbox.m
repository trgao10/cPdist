G = Mesh('off','../DATA/PNAS/meshes/s17_sas.off');
NaNInds = [2893,2980,3131];

G.DeleteVertex(2980);


figure;G.draw();
hold on;
scatter3(G.V(1,NaNInds),G.V(2,NaNInds),G.V(3,NaNInds),30,'g','filled');
G.Write('./j17_sas.off','off',[]);

