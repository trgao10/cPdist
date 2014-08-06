%Parameters

%%% **missing** TPS_disk or TPS (without partial) !!:
% Teeth : TPS, use both maxima and minima for gauss
% MT1: TPS_DISK, use both maxima and minima for gauss
% RADIUS: TPS, use local maxima for gauss points only
% base_path = '/home/grad/trgao10/MATLAB/CPsurfcompGUI/';
base_path = '/media/trgao10/Work/MATLAB/CPsurfcompGUI/';
%% dataset: choose one (uncomment)
 dataset = 'teeth';
%  dataset = 'mt1';
%  dataset = 'radius';

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% You may NOT want to modify anything below this line unless you are
% developing the codes ---------------------------------------------
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


distances_path = [base_path 'DATA/' dataset '/' 'distances'  '/'];
scripts_path = [base_path 'DATA/' dataset '/' 'scripts/'];
eo_path = [base_path 'DATA/' dataset '/' 'err_and_out/'];
mdbs_path = [base_path 'DATA/' dataset '/' 'mdbs_extrema_conf_factors_1000spreadpts'  '/'];
meshes_path = [base_path 'DATA/' dataset '/' 'meshes'  '/'];
landmarks_path = [base_path 'DATA/' dataset '/' 'landmarks/'];
cwndist_path = [base_path 'DATA/' dataset '/' 'cWn_distances/'];
data_path = [base_path 'DATA/' dataset '/'];


taxa_codes_filename = [base_path 'DATA/' dataset '/' 'taxa_codes_' dataset];

% add functions paths
path(path,[base_path 'SURFCOMP_functions']);
% path(path,[base_path 'kdtree']);

file_type = '.off';

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% PARAMETERS
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
m_or_me = 1; %m_or_me defines if we want to deform the mesh itself =1 or midedge mesh = 2
type_alg = 2; %1 - with Moser projection to area-preserving, 2- without Moser

options.num_of_spread_points = 1000; %number of density points to spread.
%(coarse representation of the surface)
%used to measure cP energy (first options.num_of_points_for_comparison)
%and Moser deformation

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%load parameters based on the current dataset >>>
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%landmark paths
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
switch dataset
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%teeth
    case 'teeth'
        landmarks_path = [landmarks_path 'landmarks_teeth'];
        
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%mt1:
    case 'mt1'
        landmarks_path = [landmarks_path 'landmarks_mt1'];
        
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % radius
    case 'radius'
        landmarks_path = [landmarks_path 'landmarks_radius'];
        
end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%method of picking the feature points for generating transofrmations
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
switch dataset
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%teeth
    case 'teeth'
        options.method_of_feature_pick = 'extrema_conf_factors'; %local extrema of conformal factors
        options.local_max_width = 5;%15; %the size of local neighborhood in searching landmarks
        options.smooth_fields = 10;%25; %number of smoothing iteration for the conformal factor (if applicable) and curvature (if applicable) field.
        options.TPS_type = 'TPS';
        
        
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % %%mt1:
    case 'mt1'
        options.method_of_feature_pick = 'mean_curvature';%'gauss_and_mean_curvature';%'gauss_curvature';%'mean_curvature';
        options.local_max_width = 5;%5;
        options.smooth_fields = 15; %10
        options.TPS_type = 'TPS_DISK';
        
    % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % radius
    case 'radius'
        options.method_of_feature_pick = 'gauss_maxima_curvature';%'gauss_and_mean_curvature';%'gauss_curvature';%'mean_curvature';
        options.local_max_width = 10;%5;
        options.smooth_fields = 10; %10
        options.TPS_type = 'TPS';
        
end


%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% For choosing seed density points and boundary layer
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

switch dataset
    
    % %  %Teeth
    case 'teeth'
        options.method_of_first_density_feature_pick = 'gauss_curvature';
        options.local_max_width_curvature = 5;
        options.smooth_curvature_fields = 10;
        options.remove_boundary_portion = 0.03; %the portion near the boundary to remove (0 - nothing)
        
    % %MT1
    case 'mt1'
        options.method_of_first_density_feature_pick = 'mean_curvature';
        options.local_max_width_curvature = 3;
        options.smooth_curvature_fields = 10;
        options.remove_boundary_portion = 0.1; %the portion near the boundary to remove (0 - nothing)
        
    % %radius
    case 'radius'
        options.method_of_first_density_feature_pick = 'gauss_maxima_curvature';
        options.local_max_width_curvature = 3;
        options.smooth_curvature_fields = 10;
        options.remove_boundary_portion = 0.05; %the portion near the boundary to remove (0 - nothing)
        
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Algorithm scanning parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

switch dataset
    
    % %%Teeth
    case 'teeth'
        options.angle_jumps = 0.05; %between 0.05-0.1 -> angle jumps [0,2pi) for the search for Mobius transformation.
        options.num_of_points_for_comparison = 50; %number of density points used to approximate the Procrustes energy
        
    % %MT1
    case 'mt1'
        options.angle_jumps = 0.1; %between 0.05-0.1 -> angle jumps [0,2pi) for the search for Mobius transformation.
        options.num_of_points_for_comparison = 50; %number of density points used to approximate the Procrustes energy
        
    %Radius
    case 'radius'
        options.angle_jumps = 0.05; %between 0.05-0.1 -> angle jumps [0,2pi) for the search for Mobius transformation.
        options.num_of_points_for_comparison = 100;%50;%200;%400; %number of density points used to approximate the Procrustes energy
        
end

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
% Cluster's jobs parameters
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

switch dataset
    
    %Teeth
    case 'teeth'
        chunk_size = 55;  %% number of comparison per job:
                            %for Teeth 161^2/X, X=number of jobs
        
    % %MT1
    case 'mt1'
        chunk_size = 15;  %%for MT1 61*61/X

        
    %Radius
    case 'radius'
        chunk_size = 10;  %%for Radius 45*45/X
        
end



