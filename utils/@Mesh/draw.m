function h = draw(G, options, FaceField)

F = G.F;
V = G.V;
if isempty(F) || isempty(V)
    disp('Empty vertices or faces!');
    return;
end
if ~iscell(F) %convert to cell array
    F=mat2cell(F',ones(1,size(F,2)),3);
end

face_size = cellfun(@numel, F); % get face size per face
minmax_face_size = [min( face_size ),max(face_size )];

if minmax_face_size(1)<3
    warning('Non polygonal faces: %s', num2str(find(fdeg<3)'));
    minmax_face_size(1) = 3;
end

h = zeros( minmax_face_size (2), 1 );

if size(F{1},1)~=1
    F = cellfun( @(x) x', F, 'UniformOutput', 0);
end

for i=minmax_face_size(1):minmax_face_size(2)
    h(i) = patch('Vertices', V', 'Faces', cell2mat(F(face_size==i)),'FaceColor',[0.5,0.5,0.5],'CDataMapping','direct');
end

h = h( h>0 );
if exist('options', 'var')==1
    set(h, options);
end

if exist('FaceField', 'var')
    faceCenter=cell2mat(cellfun(@(x)(sum(V(:,x),2)), F','UniformOutput',false))./repmat(face_size',3,1);
    hold on;
    quiver3(faceCenter(1,:),faceCenter(2,:),faceCenter(3,:),FaceField(1,:),FaceField(2,:),FaceField(3,:));
end
dim = size(V, 1);
if dim==3 && max(abs(V(3,:))) < 1e-8
    dim = 2;
end

if dim<3
    view(2);
else
    view(3);
end

grid off;
axis off;
axis tight;
axis equal; 

cameratoolbar;
cameratoolbar('SetCoordSys', 'none');

