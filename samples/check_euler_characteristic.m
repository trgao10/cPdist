samples_path = [pwd '/PNAS/'];

sampleFiles = dir(samples_path);
sampleFiles(1:2) = [];

for j=1:length(sampleFiles)
    load([samples_path sampleFiles(j).name]);
    if (G.nV-G.nE+G.nF ~= 1)
        disp([sampleFiles(j).name ' Euler Characteristic Anomaly!']);
    end
    clear G;
end
