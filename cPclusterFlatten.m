%%% preparation
clear all;
close all;
path(pathdef);
addpath(path,genpath([pwd '/utils/']));

%%% setup paths
base_path = [pwd '/'];
result_path = base_path;
data_path = '../DATA/PNAS/';

%%% clean up paths
meshes_path = [data_path 'meshes/'];
taxa_file = [data_path 'teeth_taxa_table.mat'];
samples_path = [base_path 'samples/'];
cluster_path = [base_path 'cluster/'];
scripts_path = [cluster_path 'scripts/'];
errors_path = [cluster_path 'errors/'];
outputs_path = [cluster_path 'outputs/'];

%%% load taxa codes
taxa_code = load(taxa_file);
taxa_code = taxa_code.taxa_code;
GroupSize = length(taxa_code);

%%% build folders
touch(result_path);
touch(samples_path);
touch(scripts_path);
touch(errors_path);
touch(outputs_path);

command_text = ['!rm -f ' scripts_path '*']; eval(command_text); disp(command_text);
command_text = ['!rm -f ' errors_path '*']; eval(command_text); disp(command_text);
command_text = ['!rm -f ' outputs_path '*']; eval(command_text); disp(command_text);
% command_text = ['!rm -f ' samples_path '*']; eval(command_text); disp(command_text);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(['Submitting jobs for sampling mesh files in' meshes_path '...' ]);

for k = 1:GroupSize
    job_id = k;
    
    if (exist([samples_path taxa_code{k} '_sample.mat'], 'file'))
        job_id = job_id+1;
        continue;
    end
    
    script_name = [scripts_path 'script_' num2str(job_id)];
    
    mesh_file= [meshes_path taxa_code{k} '_sas.off'];
    sample_file= [samples_path taxa_code{k} '_sample.mat'];
    
    fid = fopen(script_name,'w');
    fprintf(fid, '#!/bin/bash\n');
    fprintf(fid, '#$ -S /bin/bash\n');
    script_text = ['matlab -nodesktop -nodisplay -nojvm -nosplash -r ' ...
        ' "cd ' base_path '; ' ...
        'path(genpath(''' base_path 'utils/''), path); ' ...
        'flatten_ongrid ' ...
        mesh_file ' ' ...
        sample_file '; exit; "'];
    % system(script_text); %% grid fails on certain tasks
    fprintf(fid, '%s',script_text);
    fclose(fid);
    
    %%% qsub
    jobname = ['TCjob_' num2str(job_id)];
    serr = [errors_path 'e_job_' num2str(job_id)];
    sout = [outputs_path 'o_job_' num2str(job_id)];
    tosub = ['!qsub -N ' jobname ' -o ' sout ' -e ' serr ' ' script_name ];
    eval(tosub);
    
end
