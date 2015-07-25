function min_path = ViewRmAlongMST(name1,name2,options)
%VIEWRMALONGMST Summary of this function goes here
%   Detailed explanation goes here

if nargin<3
    options = [];
end
options.ShowTree = getoptions(options,'ShowTree','off');
options.ShowRigidMotions = getoptions(options,'ShowRigidMotions','off');
data_path = getoptions(options,'data_path','./DATA/PNAS/');

% data_path = './DATA/PNAS/';
sample_path = [data_path 'samples/'];
MST_path = [data_path 'maps_MST/'];

taxa_code = load([data_path 'teeth_taxa_table.mat']);
taxa_code = taxa_code.taxa_code;

cPdist = load([data_path 'PNAS_sym_cPdists.mat']);
cPdist = cPdist.symmetrized_all_cPdistances;
cPdist = sparse(cPdist);
cPdist = tril(cPdist, -1);

[ST, ~] = graphminspantree(cPdist, 'Method', 'Kruskal');
TAXAind1 = find(strcmpi(taxa_code, name1));
TAXAind2 = find(strcmpi(taxa_code, name2));
[~,min_path,~] = graphshortestpath(ST,TAXAind1,TAXAind2,'directed',false);
disp(['Minimal MST Path from ' name1 ' to ' name2 ': ']);
disp(taxa_code(min_path));

if strcmpi(options.ShowTree,'on')
    h = view(biograph(ST, taxa_code, 'ShowArrows', 'off', 'ShowWeights', 'on'));
    set(h.Nodes(min_path),'Color',[1 0.4 0.4])
    fowEdges = getedgesbynodeid(h,get(h.Nodes(min_path),'ID'));
    revEdges = getedgesbynodeid(h,get(h.Nodes(fliplr(min_path)),'ID'));
    edges = [fowEdges;revEdges];
    set(edges,'LineColor',[1 0 0])
    set(edges,'LineWidth',1.5)
end

if strcmpi(options.ShowRigidMotions,'on')
    PathLength = length(min_path);
    h = zeros(1,PathLength);
    if (PathLength<=5)
        NumRows = 1;
    elseif (PathLength>5 && PathLength<=10)
        NumRows = 2;
    elseif (PathLength>10 && PathLength<=15)
        NumRows = 3;
    else
        NumRows = 4;
    end
    NumCols = ceil(PathLength/NumRows);
    
    if (~isempty(findobj('Type','figure')))
        camUpVector = get(gca, 'CameraUpVector');
        camPosition = get(gca, 'CameraPosition');
        camTarget = get(gca, 'CameraTarget');
        camViewAngle = get(gca, 'CameraViewAngle');
    end
    
    figure;
    set(gcf, 'ToolBar', 'none');
    
    h(1) = subplot(NumRows, NumCols, 1);
    GM = load([sample_path name1 '_sample.mat']);
    GM = GM.GM;
    GM.draw();
    title(GM.Aux.name);
    
    chunk_size = 55;
    
    Macc = 1:max(size(GM.V));
    
    for j=2:PathLength
        locTAXAind1 = min_path(j-1);
        locTAXAind2 = min_path(j);
        GN = load([sample_path taxa_code{locTAXAind2} '_sample.mat']);
        GN = GN.GM;
        chunk_index = ceil(((locTAXAind1-1)*length(taxa_code)+locTAXAind2)/chunk_size);
        filename = [MST_path 'maps_MST_matrix_' num2str(chunk_index) '.mat'];
        maps = load(filename);
        maps = maps.maps;
        node_map = maps(locTAXAind1, locTAXAind2).cP_map_MST;
        
        R = TEETH_find_best_rigid_motion(GN.V(:,node_map)', GM.V');
        
        Macc = node_map(Macc);
        GN.V = R*GN.V;
        
        h(j) = subplot(NumRows, NumCols, j);
        GN.draw();
        title(GN.Aux.name);
        clear GM;
        GM = Mesh(GN);
        clear GN;
    end
    
    clear GM;
        
    Link = linkprop(h, {'CameraUpVector', 'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
    setappdata(gcf, 'StoreTheLink', Link);
    
    set(gca, 'CameraUpVector', [0.8469,-0.5272,-0.0696]);
    set(gca, 'CameraPosition', [0.0584,0.8255,-5.7263]);
    set(gca, 'CameraTarget', [0.0122,-0.0075,0.0173]);
    set(gca, 'CameraViewAngle', 10.5477);
    
    if (exist('camUpVector', 'var'))
        set(gca, 'CameraUpVector', camUpVector);
        set(gca, 'CameraPosition', camPosition);
        set(gca, 'CameraTarget', camTarget);
        set(gca, 'CameraViewAngle', camViewAngle);
    end
    
end


end

