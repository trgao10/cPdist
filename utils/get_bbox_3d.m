function [bbox, frame] = get_bbox_3d(pointCloud)
%GET_BBOX:    get bounding box of input point cloud
%pointCloud:  Nx3 matrix, 3 = dimension of the ambient space, N = number of
%             points in the cloud

if (size(pointCloud,1)<size(pointCloud,2))
    pointCloud = pointCloud';
end

dim = size(pointCloud,2);

if(dim~=3)
    error('This function only extract bouding box for point cloud in 3d.');
end

% [frame, ~, latent] = pca(pointCloud);
[frame, ~, latent] = princomp(pointCloud);

mark_position = zeros(dim);

for i=1:dim
    mark_position(:,i) = frame(:,i)./norm(frame(:,i));
    mark_position(:,i) = mark_position(:,i)*(max(abs(pointCloud*frame(:,i))));
end

bbox = zeros(dim, 2^dim); % dim=3
count = 1;
for i=1:2
    for j=1:2
        for k=1:2
            bbox(:,count) = mark_position*[(-1)^i, (-1)^j, (-1)^k]';
            count = count+1;
        end
    end
end

frame = frame*diag(sqrt(latent));

end
