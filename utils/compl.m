function [complX] = compl(X)
%COMPL:         complexify an 2xN matrix
%    X:         matrix of size 2xN
%    complX:    complexified X
%
%   Tingran Gao, trgao10@math.duke.edu
%   last modified: 15 Aug 2014
%

if size(X,1)>size(X,2)
    X = X';
end

complX = X(1,:)+1i*X(2,:);

end

