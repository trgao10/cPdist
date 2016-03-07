data_path = '/media/trgao10/Work/MATLAB/DATA/PNAS/';
meshes_path = [data_path 'meshes/'];

meshFiles = dir(meshes_path);
meshFiles(1:2) = [];

options.NumLandmark = 16;
for j=1:length(meshFiles)
    meshName = strtok(meshFiles(j).name,'_');
    [LandmarkInds,Landmarks] = GetLandmarks(meshName,[data_path 'landmarks_teeth.mat'],[meshes_path meshFiles(j).name],options);
    load(['./PNAS/' meshName '.mat']);
    G.Aux.Landmarks = Landmarks;
    G.Aux.LandmarkInds = LandmarkInds;
    save(['./sampels/' meshName '.mat'], 'G');
end
