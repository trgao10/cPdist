function [e_u_star, M, e2v, nume] = CORR_calculate_conjugate_harmonic_faster(F,V,mF,u,M,e2v,nume,imissing_f,reflect_mesh)
%a procedure to calculate the conjugate harmonic function u_star to u on
%the mid-edges according to polthier. This is to create conformal mapping
%of the mid edges.
% we use here th dual cot-formula

%construct the mid-edge vertex - mid-edge face ring (which has the same
%face indices as th original mesh
ring = compute_vertex_face_ring(mF);


%traverse the mesh starting for the first face
tobe_visited = zeros(nume,1);
tobe_visited(1) = 1;
tobe_len = 1;
e_u_star = NaN*ones(nume,1);
e_u_star(M(F(tobe_visited(1),1),F(tobe_visited(1),2))) = 0;%set the addditive constant

while(tobe_len>0 ) %while we havn'e finished traversing the mesh
    indf = tobe_visited(tobe_len);
    f = F(indf,:);
    tobe_len = tobe_len-1;

    %the three mid-edge vertices in this triangle
    imv1 = M(f(2),f(3));
    imv2 = M(f(3),f(1));
    imv3 = M(f(1),f(2));

%         %for debug
%     if( imv1 == 4 || imv2 == 4 || imv3 == 4)
%         disp('got to 4');
%     end
    
    
    %add to the tobe_visited only the neighboring faces which in the common edge there is no
    %data
    if (isnan(e_u_star(imv1))) %if the mid-edge vertex 1 has no data
        neigh_f=ring{imv1};
        if(neigh_f(1) == indf)
            neigh_f(1)=[];
        end
        if(~isempty(neigh_f))
            if (neigh_f(1) ~= imissing_f)%make sure we are not passing through the missing face
                tobe_len = tobe_len+1;
                tobe_visited(tobe_len) = neigh_f(1); %add this face to be visited
            end
        end
    end
    if (isnan(e_u_star(imv2))) %if the mid-edge vertex 2 has no data
        neigh_f=ring{imv2};
        if(neigh_f(1) == indf)
            neigh_f(1)=[];
        end
        if(~isempty(neigh_f))
            if (neigh_f(1) ~= imissing_f)
                tobe_len = tobe_len+1;
                tobe_visited(tobe_len) = neigh_f(1); %add this face to be visited
            end
        end
    end
    if (isnan(e_u_star(imv3))) %if the mid-edge vertex 3 has no data
        neigh_f=ring{imv3};
        if(neigh_f(1) == indf)
            neigh_f(1)=[];
        end
        if(~isempty(neigh_f))
            if (neigh_f(1) ~= imissing_f)
                tobe_len = tobe_len+1;
                tobe_visited(tobe_len) = neigh_f(1); %add this face to be visited
            end
        end
    end

    if (~isnan(e_u_star(imv1))) %if the mid-edge vertex 1 has data
        %dont change f
    elseif (~isnan(e_u_star(imv2)))
        f = [f(2) f(3) f(1)];
    elseif (~isnan(e_u_star(imv3)))
        f = [f(3) f(1) f(2)];
    else
        %not error - can happen
%         disp('Error 2341 - stopping...')
%         return
    end
        imv1 = M(f(2),f(3)); %set the midedges vertices indices again
        imv2 = M(f(3),f(1));
        imv3 = M(f(1),f(2));

%     %for debug
%     e_u_star(imv3) = 1;
%     e_u_star(imv2) = 1;
%     e_u_star(imv1) = 1;
%     if(tobe_len == 0)
%         disp('for debug');
%         break;
%     end

    
        v1 = V(f(1),:); %the regular vertices
        v2 = V(f(2),:);
        v3 = V(f(3),:);
        u1 = u(f(1)); %the discrete harmonic values st the vertices v_i
        u2 = u(f(2)); %the discrete harmonic values
        u3 = u(f(3)); %the discrete harmonic values
    
        %set the u_star at the other two vertices
%         if(reflect_mesh == 0)
            e_u_star(imv3) = CORR_set_local_conjugate_mv1_mv3(e_u_star(imv1),v1,v2,v3,u1,u2,u3);
            e_u_star(imv2) = CORR_set_local_conjugate_mv1_mv3(e_u_star(imv3),v3,v1,v2,u3,u1,u2);
%         else %should reflect the flattening
%             e_u_star(imv3) = -CORR_set_local_conjugate_mv1_mv3(-e_u_star(imv1),v1,v2,v3,u1,u2,u3);
%             e_u_star(imv2) = -CORR_set_local_conjugate_mv1_mv3(-e_u_star(imv3),v3,v1,v2,u3,u1,u2);
%         end
    
%     %for debug
%     if( imv2 == 4 || imv3 == 4)
%         disp('got to 4');
%     end

end

% % %%create edge numbering matrix: M_ji=index of undirected edge from v_i -> v_j
% % %[M e2v nume] = compute_edge_numbering(F);
% %
% % %construct the linear equations Ax=b:
% % % for each face create three equations (LS systme but should have exact
% % % solution)
% % % matrix of the system A  of size (2*#F) x #E
% % % right hand side of size (2*#F)x1
% % numf = size(F,1);
% %
% % A = sparse(2*numf,nume);
% % b = zeros(2*numf,1);
% % for k=1:numf
% %     f = F(k,:);
% %     v1 = V(f(1),:);
% %     v2 = V(f(2),:);
% %     v3 = V(f(3),:);
% %
% %     %calculate the midegde values of the harmonic function u
% %     u1 = 0.5*(u(f(2)) + u(f(3)));
% %     u2 = 0.5*(u(f(1)) + u(f(3)));
% %     u3 = 0.5*(u(f(1)) + u(f(2)));
% %
% %     %get the edges index numbers
% %     e1 = M(f(2),f(3));
% %     e2 = M(f(1),f(3));
% %     e3 = M(f(1),f(2));
% %     %calculate the angles
% %     alpha1 = myangle(v2-v1,v3-v1);
% %     alpha2 = myangle(v3-v2,v1-v2);
% %     alpha3 = myangle(v1-v3,v2-v3);
% %
% %     %write the matrix side of equations 2k-1 and 2k
% %     A(2*k-1,e3) = 1;
% %     A(2*k-1,e1) = -1;
% %
% %     A(2*k,e3) = 1;
% %     A(2*k,e2) = -1;
% %
% %     %write the right hand side
% %     b(2*k-1,1) = (cot(alpha3)*(u2-u1) + cot(alpha1)*(u2-u3));
% %     b(2*k,1) = (cot(alpha3)*(u2-u1) + cot(alpha2)*(u3-u1));
% %
% % end
% %
% % %set the constant additive
% % A(2*numf+1,1)=1;
% % b(2*numf+1) = 0.03;
% %
% % e_u_star = A\b;
% %
% % % %sample the harmonic function u on the midedges
% % % u_midp = zeros(nume,1);
% % % for k=1:nume
% % %     n12 = e2v(k,:);
% % %     n1=n12(1);
% % %     n2=n12(2);
% % %
% % %     u_midp(k) = 0.5*(u(n1)+u(n2));
% % %
% % %
% % % end
% %


