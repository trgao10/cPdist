function dVInds = DeleteIsolatedVertex(G)
%DELETEISOLATEDVERTEX Summary of this function goes here
%   Detailed explanation goes here

dVInds = find(logical((sum(G.F2V)==0)+(sum(G.F2V)==1)));
G.DeleteVertex(dVInds);

end

