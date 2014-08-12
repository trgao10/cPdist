function DeleteIsolatedVertex(G)
%DELETEISOLATEDVERTEX Summary of this function goes here
%   Detailed explanation goes here

dVInds = find(sum(G.F2V)==0);
G.DeleteVertex(dVInds);

end

