function [NA,NB]=add_dc_cr(A,B,with_bdry,mesh,~)
    [edges_n,~] = size(mesh.edges);
    is_edge_on_bdry_list = is_edge_on_bdry(mesh);
    if with_bdry
        d_idx=[];
        for i=1:size(edges_n)
            if is_edge_on_bdry_list(i)
                d_idx=[d_idx,i];
            end
        end
        A(d_idx,:)=zeros(size(A(d_idx,:)));
        A(:,d_idx)=zeros(size(A(:,d_idx)));
        B(d_idx,:)=zeros(size(B(d_idx,:)));
        B(:,d_idx)=zeros(size(B(:,d_idx)));
    else
        % delete zero cells
        d_idx=[];
        for i=1:edges_n
            if is_edge_on_bdry_list(i)
                d_idx=[d_idx,i];
            end
        end
        A(d_idx,:)=[];
        A(:,d_idx)=[];
        B(d_idx,:)=[];
        B(:,d_idx)=[];
    end
    NA = I_hull(A,A');
    NB = I_hull(B,B');

end


