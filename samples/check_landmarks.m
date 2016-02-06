Name = '32';

A = Mesh('off',['~/Work/MATLAB/DATA/Clement/meshes/' Name '.off']);
load('~/Work/MATLAB/DATA/Clement/landmarks_clement.mat');
landmarks = squeeze(PP(str2double(Name),:,:));

figure;
A.draw();
hold on;
scatter3(landmarks(:,1),landmarks(:,2),landmarks(:,3),30,'b','filled');

