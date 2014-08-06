function Write(G,filename,format,options)

options.pointCloud = getoptions(options, 'pointCloud', 0);

switch format
    case 'off'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end
        
        % header
        fprintf(fid, 'OFF\n');
        if options.pointCloud==0
            fprintf(fid, '%d %d 0\n', length(G.V), length(G.F));
        else
            fprintf(fid, '%d 0 0\n', length(G.V));
        end
        
        % write the points & faces
        fprintf(fid, '%f %f %f\n', G.V);
        if options.pointCloud==0
            fprintf(fid, '3 %d %d %d\n', G.F-1);
        end
        
        fclose(fid);
    case 'obj'
        fid = fopen(filename,'wt');
        if( fid==-1 )
            error('Can''t open the file.');
        end
        
%         if exist('options','var')==1  % open filename.mtl for writing later
%             fprintf(fid, 'mtllib %s.mtl\n', filename);
%         end
        % vertex coordinates
        fprintf(fid, 'v %f %f %f\n', G.V);
        
        % Texture coordinates
        if isfield(options, 'Texture')
            fprintf(fid, 'vt %f %f\n', options.Texture.Coordinates(1:2,:));
%             fprintf(fid, 'usemtl %s\n', filename);
            fprintf(fid, 'f %d/%d %d/%d %d/%d\n', kron(G.F',[1,1])');
        else
            fprintf(fid, 'f %d %d %d\n', G.F');
        end
        fclose(fid);
        
        
%         % MTL generation
%         if exist('options','var')==1
%             mtl_file = [filename,'.mtl'];
%             fid = fopen(mtl_file,'wt');
%             if( fid==-1 )
%                 error('Can''t open the material file.');
%             end
%             
%             Ka = [0.59 0.59 0.59];
%             Kd = [0.5 0.43 0.3];
%             Ks = [0.6 0.6 0.6];
%             d = 1;
%             Ns = 2;
%             illum = 2;
%             
%             fprintf(fid, 'newmtl %s\n', filename);
%             fprintf(fid, 'Ka  %f %f %f\n', Ka);
%             fprintf(fid, 'Kd  %f %f %f\n', Kd);
%             fprintf(fid, 'Ks  %f %f %f\n', Ks);
%             fprintf(fid, 'd  %d\n', d);
%             fprintf(fid, 'Ns  %d\n', Ns);
%             fprintf(fid, 'illum %d\n', illum);
%             fprintf(fid, 'map_Kd %s', options.Texture.Filename);
%             
%             fclose(fid);
%         end
end