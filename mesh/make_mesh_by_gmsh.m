function mesh = make_mesh_by_gmsh(a,b, h)

    global mesh_path INTERVAL_MODE

    vertices_origin = [0 0; 1 0; a b];
    vertices = [0 0; 1 0; I_mid(a) I_mid(b)];
    fid = fopen([mesh_path, 'temp.geo'],'w');
    for i = 1:3
        fprintf(fid,'Point(%d) = {%g, %g, 0, %g};\n', i, vertices(i,1), vertices(i,2), h);
    end
    fprintf(fid,'Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};\n');
    fprintf(fid,'Line Loop(1) = {1,2,3}; Plane Surface(1) = {1};\n');
    fclose(fid);
    
    %This is used for Linux/Mac OS 
    system(['bash ./mesh/create_mesh.sh ', mesh_path])
    mesh = gmshread([mesh_path, 'temp.msh']);
    %mesh = read_dolfin_mesh("/tmp/temp.xdmf")

    mesh.domain = vertices;
    if INTERVAL_MODE
        apply_exact_boundary_point_setting(mesh,vertices_origin);
    end

end
