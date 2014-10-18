samples_path = [pwd '/Clement/'];

sampleFiles = dir(samples_path);
sampleFiles(1:2) = [];

for j=1:length(sampleFiles)
    load([samples_path sampleFiles(j).name]);
    if sum(isnan(compl(G.Aux.UniformizationV)))
        disp([sampleFiles(j).name ' contains NaN in its uniformization!']);
    end
    clear G;
end

