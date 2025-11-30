function mesh = make_mesh_by_gmsh(a,b, h)
    vertices = [0 0; 1 0; a b];
    fid = fopen('temp.geo','w');
    for i = 1:3
        fprintf(fid,'Point(%d) = {%g, %g, 0, %g};\n', i, vertices(i,1), vertices(i,2), h);
    end
    fprintf(fid,'Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};\n');
    fprintf(fid,'Line Loop(1) = {1,2,3}; Plane Surface(1) = {1};\n');
    fclose(fid);

    system('/opt/homebrew/bin/gmsh temp.geo -2 -format msh2 -o temp.msh');

    mesh = gmshread('temp.msh');
    mesh.domain = vertices;
end