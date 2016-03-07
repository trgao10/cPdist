clear vars;
path(pathdef);
addpath(path,genpath([pwd '/../utils/']));

% data_path = '/media/trgao10/Work/MATLAB/ArchivedData/131astragali_52species/';
% data_path = '/media/trgao10/Work/MATLAB/DATA/PNAS/meshes/';
data_path = '/media/trgao10/Work/MATLAB/DATA/Clement/meshes/';
[ds.names, suffix] = getFileNames(data_path);

for j=1:length(ds.names)
    progressbar(j, length(ds.names), 20);
    [Ver, F] = read_off([data_path ds.names{j} suffix]);
%     [Ver, F] = read_ply([data_path ds.names{j} suffix]);
    Ver = Ver - repmat(mean(Ver), size(Ver, 1), 1);
    [U,~,~] = svd(Ver);
    if det(U) < 0
        disp([ds.names{j} ' flipped']);
    end
end


% spreadsheet_path = [data_path 'ClassificationTable.xlsx'];
% sample_path = '/home/trgao10/Downloads/cPdistMinimal/PNAS/';

% taxa_file = [data_path 'teeth_taxa_table.mat'];
% taxa_code = load(taxa_file);
% taxa_code = taxa_code.taxa_code;

% TB = readtable([data_path 'ClassificationTable.xlsx']);
% Name = TB.Name;
% Chirality = TB.Chirality;

% idxL = find(strcmpi(Chirality, 'Left'));
% idxR = find(strcmpi(Chirality, 'Right'));
% 
% for j=1:length(idxL)
%     idxL(j) = find(strcmpi(taxa_code, Name{idxL(j)}));
% end
% for j=1:length(idxR)
%     idxR(j) = find(strcmpi(taxa_code, Name{idxR(j)}));
% end
% 
% Gs = cell(1, length(taxa_code));
% for j=1:length(taxa_code)
%     load([sample_path taxa_code{j} '.mat']);
%     Gs{j} = G;
% end
% 
% for j=1:length(taxa_code)
%     progressbar(j, length(taxa_code), 20);
%     [UX,DX,VX] = svd(Gs{j}.V);
%     if det(UX) < 0
%         disp([Gs{j}.Aux.name ' reverted']);
%     end
% end

% disp('start searching left teeth');
% for j=1:length(idxL)
%     [UX,DX,VX] = svd(Gs{idxL(j)}.V);
%     for k = (j+1):length(idxL)
%         [UY,DY,VY] = svd(Gs{idxL(k)}.V);
%         if det(UX*UY') < 0
%             disp([Gs{idxL(j)}.Aux.name ' and ' Gs{idxL(k)}.Aux.name ' not consistent']);
%         end
%     end
% end
% 
% disp('start searching right teeth');
% for j=1:length(idxR)
%     [UX,DX,VX] = svd(Gs{idxR(j)}.V);
%     for k = (j+1):length(idxR)
%         [UY,DY,VY] = svd(Gs{idxR(k)}.V);
%         if det(UX*UY') < 0
%             disp([Gs{idxr(j)}.Aux.name ' and ' Gs{idxr(k)}.Aux.name ' not consistent']);
%         end
%     end
% end


% for j=1:length(idxL)
%     [UX,DX,VX] = svd(Gs{idxL(j)}.V);
%     for k = 1:length(idxR)
%         [UY,DY,VY] = svd(Gs{idxR(k)}.V);
%         
%         if det(UX*UY') > 0
%             disp([Gs{idxL(j)}.Aux.name ' and ' Gs{idxR(k)}.Aux.name ' not consistent']);
%         end
%     end
% end
