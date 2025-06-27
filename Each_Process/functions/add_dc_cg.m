function [NA,NB]=add_dc_cg(A,B,with_bdry,mesh,~)
    [nodes_n,~]=size(mesh.nodes);
    is_node_on_bdry_list = is_node_on_bdry(mesh);

    if with_bdry
        d_idx=[];
        for i=1:nodes_n
            if is_node_on_bdry_list(i)
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
        for i=1:nodes_n
            if is_node_on_bdry_list(i)
                d_idx=[d_idx,i];
            end
        end
        A(d_idx,:)=[];
        A(:,d_idx)=[];
        B(d_idx,:)=[];
        B(:,d_idx)=[];
    end
    NA = (A+A')/2;
    NB = (B+B')/2;
end


