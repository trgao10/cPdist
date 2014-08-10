function [LandmarkInds,Landmarks] = GetLandmarks(G,LandmarksPath)
%GETLANDMARKS Summary of this function goes here
%   Detailed explanation goes here

NumLandmark = 16;
LandmarkFile = load(LandmarksPath);
rawLandmarks = LandmarkFile.PP(strcmpi(LandmarkFile.names,G.Aux.name),1:NumLandmark,:);
Landmarks = zeros(size(rawLandmarks,2),3);
for k=1:size(rawLandmarks,2)
    Landmarks(k,:) = [rawLandmarks(1,k,1), rawLandmarks(1,k,2), rawLandmarks(1,k,3)];
end
Landmarks = Landmarks-repmat(G.Aux.Center',NumLandmark,1);
tree = KDTreeSearcher(G.V');
LandmarkInds = tree.knnsearch(Landmarks);

end

