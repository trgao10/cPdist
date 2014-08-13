%%% preparation
clear all;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
rslts_path = [base_path 'rslts/'];

%%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
chunk_size = 55;

%%% read rslt matrices and separate distance and landmarkMSE's
cPdistMatrix = zeros(GroupSize,GroupSize);
cPmapsMatrix = cell(GroupSize,GroupSize);
invcPmapsMatrix = cell(GroupSize,GroupSize);
lmkMSEMatrix = zeros(GroupSize,GroupSize);

cnt = 0;
job_id = 0;
for k1=1:GroupSize
    progressbar(k1,GroupSize,20);
    for k2=1:GroupSize
        if mod(cnt,chunk_size)==0
            job_id = job_id+1;
            load([rslts_path 'rslt_mat_' num2str(job_id)]);
        end
        cPdistMatrix(k1,k2) = cPrslt{k1,k2}.cPdist;
        cPmapsMatrix{k1,k2} = cPrslt{k1,k2}.cPmap;
        invcPmapsMatrix{k1,k2} = cPrslt{k1,k2}.invcPmap;
        lmkMSEMatrix(k1,k2) = cPrslt{k1,k2}.lkMSE;
        
        cnt = cnt+1;
    end
end

%%% symmetrize
for j=1:GroupSize
    progressbar(j,GroupSize,20);
    for k=1:GroupSize
        if cPdistMatrix(j,k)<cPdistMatrix(k,j)
            lmkMSEMatrix(k,j) = lmkMSEMatrix(j,k);
            cPmapsMatrix{k,j} = invcPmapsMatrix(j,k);
        else
            lmkMSEMatrix(j,k) = lmkMSEMatrix(k,j);
            cPmapsMatrix{j,k} = invcPmapsMatrix(k,j);
        end
    end
end
cPdistMatrix = min(cPdistMatrix,cPdistMatrix');

%%% visualize distance and landmarkMSE matrices
figure;
imagesc(cPdistMatrix./max(cPdistMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

figure;
imagesc(lmkMSEMatrix./max(lmkMSEMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

%%% save results
save('cPdistMatrix','cPdistMatrix');
save('lmkMSEMatrix','lmkMSEMatrix');
save('cPmapsMatrix','cPmapsMatrix');

