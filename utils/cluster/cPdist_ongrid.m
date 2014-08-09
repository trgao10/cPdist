function [best_Dist flat1 flat2 map_12 best_corr] = TEETH_compare_disk_surfaces_cluster( ...
                                mdb_file_1, mdb_file_2, ...
                                job_id)
% load parameters 
parameters;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the taxa

load(taxa_codes_filename);
taxa = taxa_code;

% %for doug :
% load('taxa_code_TEETH');
%  taxa=textdata; 

% %for mt1 :
%taxa_codes = load('mt1_taxa_codes');
%taxa = taxa_codes.taxa_codes; %this is the way the matrix is indexed - this is the master taxa list

%  %for biren :
%  taxa_codes = load('taxa_code_BIREN');
%  taxa = taxa_codes.taxa_code; %this is the way the matrix is indexed - this is the master taxa list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compare the surfaces
[best_Dist all_dist flat1 flat2 map_12 best_corr best_a best_tet best_ref] = ...
    TEETH_compare_disk_surfaces_kdtreesearcher([mdbs_path mdb_file_1], ...
                                [mdbs_path mdb_file_2], type_alg ,m_or_me );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load the matrix of this jobid (to fill the new entry)
dist_mat_filename = [distances_path 'dist_mat_' num2str(job_id)];
load(dist_mat_filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the indices of the currently compared teeth (by matching to taxa)

% %for doug:
% tname1= lower(filename_mdb1(1:(end-8)));tname1(1:7)=[];
% tname2 = lower(filename_mdb2(1:(end-8)));tname2(1:7)=[];

% %for biren:
% filename_mdb1(1:13)=[];
% filename_mdb2(1:13)=[];

% %for mt1:
% mdb_file_1(1:11)=[];
% mdb_file_2(1:11)=[];


tname1= strtok(lower(mdb_file_1),['_' '.']) ;
tname2 = strtok(lower(mdb_file_2),['_' '.']); 

TAXAind1 = find(strcmp(lower(taxa),lower(tname1)))
TAXAind2 = find(strcmp(lower(taxa),lower(tname2) ))

%for debug
taxa;
tname1
tname2
if(isempty(TAXAind1) || isempty(TAXAind2))
    disp('error: didnt find a taxa!');
end

%output results in the dist matrix and save it back
dist(TAXAind1,TAXAind2) = best_Dist;
save(dist_mat_filename,'dist');

