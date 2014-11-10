DistMatrix = load('./results/Teeth/cPdist/cPDistMatrix.mat');
DistMatrix = DistMatrix.cPDistMatrix;
DistMatrix = DistMatrix - diag(diag(DistMatrix));

Y = mdscale(DistMatrix,2);
plot(Y(:,1),Y(:,2),'bo');
axis equal;