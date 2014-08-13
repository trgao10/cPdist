function min_path = FindMSTPath(TAXAind1,TAXAind2,cPdistMatrix)
%FINDMSTAPTH Summary of this function goes here
%   Detailed explanation goes here

cPdist = sparse(cPdistMatrix);
cPdist = tril(cPdist, -1);

[ST, ~] = graphminspantree(cPdist, 'Method', 'Kruskal');
[~,min_path,~] = graphshortestpath(ST,TAXAind1,TAXAind2,'directed',false);

end

