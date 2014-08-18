function [ST] = ConstructGraph(DistMatrix,Type)
%CONSTRUCTGRAPH Summary of this function goes here
%   Detailed explanation goes here

if strcmpi(Type,'MST')
    DistMatrix = sparse(DistMatrix);
    DistMatrix = tril(DistMatrix, -1);
    [ST, ~] = graphminspantree(DistMatrix, 'Method', 'Kruskal');
elseif strcmpi(Type,'SPT')
    DistMatrix = sparse(DistMatrix);
    GroupSize = size(DistMatrix,1);
    for j=1:GroupSize
        [dist, path, pred] = graphshortestpath(DistMatrix, j);
    end
elseif strcmpi(Type,'Viterbi')
else
    error('Unidefined Graph Type!');
end

end

