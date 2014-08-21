%%% preparation
clear all;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% set ImprType
ImprType = 'MST';

%%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
result_path = '/xtmp/MATLAB/cPMST/'; 
rslts_path = [result_path 'rslts/'];
TextureCoords1Matrix_path = [result_path 'TextureCoords1Matrix/'];
TextureCoords2Matrix_path = [result_path 'TextureCoords2Matrix/'];

%%% check if paths exist
touch(TextureCoords1Matrix_path);
touch(TextureCoords2Matrix_path);

%%% load taxa codes
taxa_file = [data_path 'teeth_taxa_table.mat'];
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);
chunk_size = 55;

%%% read rslt matrices and separate distance and landmarkMSE's
ImprDistMatrix = zeros(GroupSize,GroupSize);
ImprMapsMatrix = cell(GroupSize,GroupSize);
invImprMapsMatrix = cell(GroupSize,GroupSize);
lmkMSEMatrix = zeros(GroupSize,GroupSize);
TextureCoords1Matrix = cell(GroupSize,GroupSize);
TextureCoords2Matrix = cell(GroupSize,GroupSize);

cnt = 0;
job_id = 0;
for k1=1:GroupSize
    progressbar(k1,GroupSize,20);
    for k2=1:GroupSize
        if mod(cnt,chunk_size)==0
            job_id = job_id+1;
            load([rslts_path 'rslt_mat_' num2str(job_id)]);
        end
        ImprDistMatrix(k1,k2) = Imprrslt{k1,k2}.ImprDist;
        ImprMapsMatrix{k1,k2} = Imprrslt{k1,k2}.ImprMap;
        invImprMapsMatrix{k1,k2} = Imprrslt{k1,k2}.invImprMap;
        lmkMSEMatrix(k1,k2) = Imprrslt{k1,k2}.lkMSE;
        TextureCoords1Matrix{k1,k2} = Imprrslt{k1,k2}.TextureCoords1;
        TextureCoords2Matrix{k1,k2} = Imprrslt{k1,k2}.TextureCoords2;
        
        cnt = cnt+1;
    end
end

%%% symmetrize
for j=1:GroupSize
    progressbar(j,GroupSize,20);
    for k=1:GroupSize
        if ImprDistMatrix(j,k)<ImprDistMatrix(k,j)
            lmkMSEMatrix(k,j) = lmkMSEMatrix(j,k);
            ImprMapsMatrix{k,j} = invImprMapsMatrix{j,k};
            TextureCoords1Matrix{k,j} = TextureCoords2Matrix{j,k};
            TextureCoords2Matrix{k,j} = TextureCoords1Matrix{j,k};
        else
            lmkMSEMatrix(j,k) = lmkMSEMatrix(k,j);
            ImprMapsMatrix{j,k} = invImprMapsMatrix{k,j};
            TextureCoords1Matrix{j,k} = TextureCoords2Matrix{k,j};
            TextureCoords2Matrix{j,k} = TextureCoords1Matrix{k,j};
        end
    end
end
ImprDistMatrix = min(ImprDistMatrix,ImprDistMatrix');

%%% visualize distance and landmarkMSE matrices
figure;
imagesc(ImprDistMatrix./max(ImprDistMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

figure;
imagesc(lmkMSEMatrix./max(lmkMSEMatrix(:))*64);
axis equal;
axis([1,GroupSize,1,GroupSize]);

%%% save results
save([result_path '/cP' ImprType 'DistMatrix'],'ImprDistMatrix');
save([result_path '/cP' ImprType 'lmkMSEMatrix'],'lmkMSEMatrix');
save([result_path '/cP' ImprType 'MapsMatrix'],'ImprMapsMatrix');

