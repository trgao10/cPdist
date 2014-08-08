function [pmV ] = CORR_map_mesh_to_plane_nonconforming(V,F,mF,seed_face,M,E2V,numE, reflect_mesh)
% map the mesh V,F to plane using the non-conforming method of polthier
% use the seed face to cut the mesh open
% output: the midpoint mesh and its embedding to plane

oF = F; %rememeber the face list for later

%cut out the seed face
ioutf = seed_face;
outf = F(ioutf,:);
outf=sort(outf);
F(ioutf,:)=[];

[L] = CORR_compute_laplacian_tension(V,F);
% %if L has row with negative entries change to combinatorial laplacian
% cnt=0;
% for k=1:size(L,1)
%    tind = find(L(k,:)<0);
%    if (~isempty(tind))%there are negative elements
%        tind = find(L(k,:));
%        m = length(tind);
%        L(k,k)=-m;
%        L(k,tind) = 1;
%        cnt=cnt+1;       
%    end
% end
% str = ['corrected ' num2str(cnt) ' row out of ' num2str(size(L,1))];
% disp(str);
% %end correcting laplcian matrix

%[L] = compute_mvc_laplacian(V,F);
L1=L;
outf=sort(outf);
L1(outf(1),:)=[];
L1(outf(2)-1,:)=[];

% %% consider pole to be centroid of the cut-out triangle
% % assign values to this triangle based on the expected pole at that point
% L1(outf(3)-2,:)=[];
% L1rows = size(L1,1);

% v1 = V(outf(1),:);
% v2 = V(outf(2),:);
% v3 = V(outf(3),:);
% c = (v1+v2+v3)/3;
% c_v1 = v1-c;
% c_v2 = v2-c;
% c_v3 = v3-c;
% L1(L1rows+1,outf(1))=1;
% L1(L1rows+2,outf(2))=1;
% L1(L1rows+3,outf(3))=1;
% 
% b=zeros(L1rows,1);
% % % take equilateral
% b(L1rows+1) = c_v1(1) / dot(c_v1,c_v1);
% b(L1rows+2) = c_v2(1) / dot(c_v2,c_v2);
% b(L1rows+3) = c_v3(1) / dot(c_v3,c_v3);
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Option 2: just place two vertices
L1rows = size(L1,1);

L1(L1rows+1,outf(1))=1;
L1(L1rows+2,outf(2))=1;

b=zeros(L1rows,1);
% % take equilateral
% if(reflect_mesh == 0)
    b(L1rows+1) = -1;
    b(L1rows+2) = 1;
% else
%     b(L1rows+1) = 1;
%     b(L1rows+2) = -1;
% end
    






u = L1\b;

% %with linear system
% [e_u_star] = CORR_calculate_conjugate_harmonic(F,V,u,M,E2V,numE);
%withOUT linear system
imissing_f = seed_face;
[e_u_star] = CORR_calculate_conjugate_harmonic_faster(oF,V,mF,u,M,E2V,numE,imissing_f,reflect_mesh);

%[mV mF] = create_midpoint_mesh(V,F);
[mu] = CORR_get_midpoint_values_of_function(oF,u, M, E2V, numE);

if(reflect_mesh==0)
    pmV = [mu e_u_star]; 
else
    pmV = [mu -e_u_star]; 
end

