function mesh = make_mesh_by_gmsh(a,b, h)

    global mesh_path

    vertices = I_intval([0 0; 1 0; a b]);
    fvertices = [0 0; 1 0; I_mid(a) I_mid(b)];
    fid = fopen([mesh_path, 'temp.geo'],'w');

    
    for i = 1:3
        fprintf(fid,'Point(%d) = {%g, %g, 0, %g};\n', i, fvertices(i,1), fvertices(i,2), I_mid(I_intval(h)));
    end
    fprintf(fid,'Line(1)={1,2}; Line(2)={2,3}; Line(3)={3,1};\n');
    fprintf(fid,'Line Loop(1) = {1,2,3}; Plane Surface(1) = {1};\n');
    fclose(fid);
    
    %This is used for Linux 
    % system(['bash ./mesh/create_mesh.sh ', mesh_path])
    [status, cmdout] = system(['bash ./mesh/create_mesh.sh ', mesh_path]);
    
    mesh = gmshread([mesh_path, 'temp.msh']);
    mesh.domain = vertices;

    mesh = apply_exact_boundary_point_setting(mesh, vertices);
    
end
