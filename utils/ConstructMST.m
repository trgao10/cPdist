function [ST] = ConstructMST(DistMatrix,TaxaCode)
%CONSTRUCTMST Summary of this function goes here
%   Detailed explanation goes here

DistMatrix = sparse(DistMatrix);
DistMatrix = tril(DistMatrix, -1);

[ST, ~] = graphminspantree(DistMatrix, 'Method', 'Kruskal');

h = view(biograph(ST, TaxaCode, 'ShowArrows', 'off', 'ShowWeights', 'on'));

end

