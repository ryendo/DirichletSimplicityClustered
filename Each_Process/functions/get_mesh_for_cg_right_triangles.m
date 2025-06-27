function mesh = get_mesh_for_cg_right_triangles(varargin)
        
    % the edge OA must be parallel to the x-axis
    tri_node= varargin{1};
    N = varargin{2};

    r=sqrt(tri_node(3,1)^2+tri_node(3,2)^2);
    th=atan(tri_node(3,2)/tri_node(3,1));

    % scatter(mid(tri_node(:,1)),mid(tri_node(:,2)));

    p2=aff([tri_node(1,1),tri_node(1,2)],r,th);
    p3=aff([tri_node(2,1),tri_node(2,2)],r,th);
    p1=aff([tri_node(3,1),tri_node(3,2)],r,th);
    
    
    x1 = p1(1); y1 = p1(2);
    x2 = p2(1); y2 = p2(2);
    x3 = p3(1); y3 = p3(2);
    
%   make nodes
    nodes=intval(zeros((N+1)^2,2));
%   left side
    dx=abs(x3-x1)/N;
    dy=abs(y3-y1)/N;
    
    for row=1:N+1
        for col=1:row
            x=x3+dx*(col-1)-dx*(row-1);
            y=y3-dy*(row-1);
            idx=row^2-2*row+2+(col-1);
            nodes(idx,:)=[x,y];
        end
    end
    
%   right side
    dx=abs(x2-x3)/N;
    dy=abs(y3-y2)/N;
    
    for row=2:N+1
        for col=2:row
            x=x3+dx*(col-1);
            y=y3-dy*(row-1);
            idx=row^2-row+1+(col-1);
            nodes(idx,:)=[x,y];
        end
    end

    [N_nodes,~]=size(nodes);
    for idx=1:N_nodes
        nodes(idx,:)=invaff(nodes(idx,:),r,th);
    end

    
    
%   set element number
    elements=[];
%   left side
    idx=1;
    elements(idx,:)=[1,2,3];
    idx=idx+1;
    for row=2:N
        a=row^2-2*row+2;
        b=(row+1)^2-2*(row+1)+2;
        c=b+1;
        elements(idx,:)=[a,b,c]; 
        idx=idx+1;
        for col=2:row
            a=row^2-2*row+2+(col-2);
            b=a+1;
            c=(row+1)^2-2*(row+1)+1+col;
            d=c+1;
            elements(idx,:)=[b,a,c];
            idx=idx+1;
            elements(idx,:)=[b,c,d];
            idx=idx+1;
        end
    end
%   right side
    elements(idx,:)=[1,3,4];
    idx=idx+1;
    for row=2:N
        for col=1:row-1
            a=row^2-row+1+(col-1);
            b=a+1;
            c=(row+1)^2-row+col-1;
            d=c+1;
            elements(idx,:)=[a,c,d];
            idx=idx+1;
            elements(idx,:)=[a,d,b];
            idx=idx+1;
        end
        elements(idx,:)=[row^2,(row+1)^2-1,(row+1)^2];
        idx=idx+1;
    end

%   make edge
    edges=[];
    idx=1;
%   left side
    edges(idx,:)=elements(1,[1,2]); idx=idx+1; edges(idx,:)=elements(1,[2,3]); idx=idx+1; edges(idx,:)=elements(1,[3,1]); idx=idx+1;
    for row=2:N
        for col=1:row
            el=row^2-2*row+2+2*(col-1);
            edges(idx,:)=elements(el,[1,2]); idx=idx+1; edges(idx,:)=elements(el,[2,3]); idx=idx+1; edges(idx,:)=elements(el,[3,1]); idx=idx+1;
        end
    end
%   right side
    edges(idx,:)=elements(N^2+1,[2,3]); idx=idx+1; edges(idx,:)=elements(N^2+1,[3,1]); idx=idx+1;
    for row=2:N
        el=N^2+row^2-2*row+2;
        edges(idx,:)=elements(el,[2,3]); idx=idx+1; edges(idx,:)=elements(el,[3,1]); idx=idx+1;
        for col=2:row
            el=N^2+row^2-2*row+2+2*(col-1);
            edges(idx,:)=elements(el,[1,2]); idx=idx+1; edges(idx,:)=elements(el,[2,3]); idx=idx+1; edges(idx,:)=elements(el,[3,1]); idx=idx+1;
        end
    end
    
% make domain
    domain=[invaff(p1,r,th);invaff(p3,r,th);invaff(p2,r,th)];
    mesh=struct('nodes',nodes,'edges',edges,'elements',elements,'domain',domain);
end

function tx=aff(x,r,th)
    th=intval(pi)-th;
    tx=1/r*[cos(th),-sin(th);sin(th),cos(th)]*x'+[1;0];
    tx=tx';
end

function x=invaff(tx,r,th)
    th=intval(pi)-th;
    x=r*[cos(-th),-sin(-th);sin(-th),cos(-th)]*(tx'-[1;0]);
    x=x';
end