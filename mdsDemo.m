DistMatrix = load('./results/Clement/cPComposedLAST/FeatureFixOn/cPComposedLASTDistMatrix.mat');
DistMatrix = DistMatrix.ImprDistMatrix;
DistMatrix = DistMatrix - diag(diag(DistMatrix));

Y = mdscale(DistMatrix,2);
plot(Y(:,1),Y(:,2),'bo');
axis equal;