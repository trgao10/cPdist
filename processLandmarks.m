meshPath = '../DATA/Clement/meshes/';
landmarkPath = '../DATA/Clement/landmarks/';
NumLandmarks = 7;

dir_struct = dir(meshPath);
names = cell(length(dir_struct)-2,1);
PP = zeros(length(dir_struct)-2,NumLandmarks,3);

for j=3:length(dir_struct)
    if (dir_struct(j).isdir~=1)
        meshName = strtok(dir_struct(j).name, '.');
        names{j-2} = meshName;
    end
    fid = fopen([landmarkPath num2str(str2double(meshName)) '.landmarkAscii']);
    if (fid==-1)
        error(['Can''t open file' landmarkPath num2str(str2double(meshName)) '.landmarkAscii']);
    end
    count = 0;
    tline = fgetl(fid);
    while ~feof(fid);
        if (numel(tline)>0 && ismember(tline(1),'0123456789'))
            count = count+1;
            PP(j-2,count,:) = sscanf(tline, '%f %f %f');
        end
        tline = fgetl(fid);
    end
    fclose(fid);
end

save([landmarkPath 'landmarks_clement.mat'], 'names', 'PP');
