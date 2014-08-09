%script: compare all pairs of surfaces 
parameters;
disp('++++++++++++++++++++++++++++++++++++++++++++++++++');
disp(['Submitting jobs for comparing all pairs in dataset ' dataset ' ... ']);
disp(' ');


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
type = '.mat';

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% distances files paths

%%clean these paths
command_text = ['!rm -f ' distances_path '*']; eval(command_text);
disp(command_text);
command_text = ['!rm -f ' scripts_path '*']; eval(command_text);
disp(command_text);
command_text = ['!rm -f ' eo_path '*']; eval(command_text);
disp(command_text);

dir_struct = dir(mdbs_path);
[mdbs_filenames] = CORR_find_all_meshes_in_dir(dir_struct,type);


% load taxa codes
taxa_codes = load(taxa_codes_filename);
taxa = taxa_codes.taxa_code;
ntaxa=length(taxa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%for every pairs of mdbs compute the distance between them (w/wo
%%reflection)
n = length(mdbs_filenames);
cnt=0;
job_id=0; %this is the job_id
for k1=1:n
    for k2 =1:n
        if( mod(cnt,chunk_size) == 0)  %save the old script and create new script
            
            if (job_id > 0 ) %not the first time
                %%% close the script file (except the last one, see below)
                fprintf(fid, 'exit; "\n', script_text);
                fclose(fid);
                
                %%% qsub
                jobname = ['TCjob_' num2str(job_id)];
                serr = [eo_path 'e_job_' num2str(job_id)];
                sout = [eo_path 'o_job_' num2str(job_id)];
                tosub = ['!qsub -q all.q -N ' jobname ' -o ' sout ' -e ' serr ' ' script_one_dist_name ];
                eval(tosub);
            end
            
            job_id = job_id+1;
            script_one_dist_name = [scripts_path ...
                'script_' num2str(job_id)];
            
            %%% open the next (first?) script file
            fid = fopen(script_one_dist_name,'w');
            fprintf(fid, '#!/bin/bash\n');
            fprintf(fid, '#$ -S /bin/bash\n');
            script_text = ['matlab -nodesktop -nodisplay -nojvm -nosplash -r '...
                '" cd ' base_path ' ; parameters; ' ];
            fprintf(fid, '%s ',script_text);
            
            %%% create new matrix
            dist=zeros(ntaxa,ntaxa);
            save([distances_path 'dist_mat_' num2str(job_id)],'dist');
            %disp(['saved: ' distancespath 'dist_mat_' num2str(job_id)]);
            
        end
        filename1 = mdbs_filenames(k1).name;
        filename2 = mdbs_filenames(k2).name;
        
        %%% prepare the script
        
        %%% using comparison with tps deformation
        script_text = [' TEETH_compare_disk_surfaces_cluster ' ...
            filename1 ' ' ...
            filename2  ' ' ...
            num2str(job_id) '; ' ];
        fprintf(fid, '%s ',script_text);
        
        cnt=cnt+1;
        
    end
    
end

%%% close the last script file
fprintf(fid, 'exit; "\n',script_text);
fclose(fid);

%%% qsub
jobname = ['TCjob_' num2str(job_id)];
serr = [eo_path 'e_job_' num2str(job_id)];
sout = [eo_path 'o_job_' num2str(job_id)];
tosub = ['!qsub -q all.q -N ' jobname ' -o ' sout ' -e ' serr ' ' script_one_dist_name ];
eval(tosub);


