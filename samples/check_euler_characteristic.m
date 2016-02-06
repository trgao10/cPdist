samples_path = [pwd '/HDM/'];
% samples_path = '/media/trgao10/Work/MATLAB/DATA/HDM/samples/';

sampleFiles = dir(samples_path);
sampleFiles(1:2) = [];

for j=1:length(sampleFiles)
    load([samples_path sampleFiles(j).name]);
    if (G.nV-G.nE+G.nF ~= 1)
        disp([sampleFiles(j).name ' Euler Characteristic Anomaly!']);
        G.DeleteIsolatedVertex();
        disp(['after removing isolated vertex, Euler Characteristics = ' num2str(G.nV-G.nE+G.nF)]);
    end
    clear G;
end
