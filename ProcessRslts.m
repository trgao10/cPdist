%%% preparation
clear all;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% setup paths
base_path = [pwd '/'];
data_path = '../DATA/PNAS/';
rslts_path = '/media/trgao10/Work/MATLAB/cPdist/rslts/';
% rslts_path = [base_path 'rslts/'];
TextureCoords1Matrix_path = '/media/trgao10/Work/MATLAB/cPdist/TextureCoords1Matrix/';
TextureCoords2Matrix_path = '/media/trgao10/Work/MATLAB/cPdist/TextureCoords2Matrix/';

%%% check if paths exist
touch(TextureCoords1Matrix_path);
touch(TextureCoords2Matrix_path);

%%% clean up texture coordinates matrices
command_text = ['!rm -f ' TextureCoords1Matrix_path '*'];
eval(command_text); disp(command_text);
command_text = ['!rm -f ' TextureCoords2Matrix_path '*'];
eval(command_text); disp(command_text);

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
tmpTextureCoords1Matrix = cell(GroupSize,GroupSize);
tmpTextureCoords2Matrix = cell(GroupSize,GroupSize);

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
        tmpTextureCoords1Matrix{k1,k2} = cPrslt{k1,k2}.TextureCoords1;
        tmpTextureCoords2Matrix{k1,k2} = cPrslt{k1,k2}.TextureCoords2;
        
        cnt = cnt+1;
    end
end

%%% symmetrize
cnt = 0;
job_id = 0;
for j=1:GroupSize
    progressbar(j,GroupSize,20);
    for k=1:GroupSize
        if mod(cnt,chunk_size)==0
            if cnt>0
                save([TextureCoords1Matrix_path 'TextureCoords1Matrix_' num2str(job_id) '.mat'],'TextureCoords1Matrix');
                save([TextureCoords2Matrix_path 'TextureCoords2Matrix_' num2str(job_id) '.mat'],'TextureCoords2Matrix');
                clear TextureCoords1Matrix TextureCoords2Matrix
            end
            job_id = job_id+1;
            TextureCoords1Matrix = cell(GroupSize,GroupSize);
            TextureCoords2Matrix = cell(GroupSize,GroupSize);
        end
        if cPdistMatrix(j,k)<cPdistMatrix(k,j)
            lmkMSEMatrix(k,j) = lmkMSEMatrix(j,k);
            cPmapsMatrix{k,j} = invcPmapsMatrix{j,k};
            TextureCoords1Matrix{j,k} = tmpTextureCoords1Matrix{j,k};
            TextureCoords2Matrix{j,k} = tmpTextureCoords2Matrix{j,k};
%             TextureCoords1Matrix{k,j} = tmpTextureCoords2Matrix{j,k};
%             TextureCoords2Matrix{k,j} = tmpTextureCoords1Matrix{j,k};
        else
            lmkMSEMatrix(j,k) = lmkMSEMatrix(k,j);
            cPmapsMatrix{j,k} = invcPmapsMatrix{k,j};
%             TextureCoords1Matrix{k,j} = tmpTextureCoords1Matrix{k,j};
%             TextureCoords2Matrix{k,j} = tmpTextureCoords2Matrix{k,j};
            TextureCoords1Matrix{j,k} = tmpTextureCoords2Matrix{k,j};
            TextureCoords2Matrix{j,k} = tmpTextureCoords1Matrix{k,j};
        end
        cnt = cnt+1;
    end
end
if mod(cnt,chunk_size)~=0
    save([TextureCoords1Matrix_path 'TextureCoords1Matrix_' num2str(job_id) '.mat'],'TextureCoords1Matrix');
    save([TextureCoords2Matrix_path 'TextureCoords2Matrix_' num2str(job_id) '.mat'],'TextureCoords2Matrix');
    clear TextureCoords1Matrix TextureCoords2Matrix
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

