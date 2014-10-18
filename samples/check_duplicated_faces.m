samples_path = [pwd '/Clement/'];

sampleFiles = dir(samples_path);
sampleFiles(1:2) = [];

for j=1:length(sampleFiles)
    load([samples_path sampleFiles(j).name]);
    [F,Inds] = unique(sort(G.F',2),'rows','first');
    if ~isempty(setdiff(1:G.nF,Inds))
        disp([sampleFiles(j).name ' contains duplicate faces!']);
    end
    clear G;
end

