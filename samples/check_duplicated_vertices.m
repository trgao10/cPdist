samples_path = [pwd '/HDM/'];

sampleFiles = dir(samples_path);
sampleFiles(1:2) = [];

for j=1:length(sampleFiles)
    load([samples_path sampleFiles(j).name]);
    [F,Inds] = unique(sort(G.V',2),'rows','first');
    if ~isempty(setdiff(1:G.nV,Inds))
        disp([sampleFiles(j).name ' contains duplicate vertices!']);
    end
    clear G;
end

