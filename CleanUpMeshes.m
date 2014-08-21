meshesPath = '../DATA/PNAS/meshes/';

meshFiles = dir(meshesPath);
meshFiles(1:2) = [];

options.display = 'on';
options.exclude_boundary = 1;

for j=1:length(meshFiles)
    progressbar(j,length(meshFiles),20);
    if ~strcmpi(meshFiles(j).name,'V13_sas.off')
        continue;
    end
    keyboard
    G = Mesh('off',[meshesPath meshFiles(j).name]);
    dVInds = G.DeleteIsolatedVertex(options);
    if ~isempty(dVInds);
        disp([meshFiles(j).name ' contains non-boundary isolated vertex!']);
        G.Write([newMeshesPath meshFiles(j).name],'off',[]);
    end
    clear G;
end

