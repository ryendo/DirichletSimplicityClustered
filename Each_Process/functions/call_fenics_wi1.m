function [p_h_matrix, g_i_h_matrix] = call_fenics_wi1(mesh)
    % Extract mesh data
    nodes = mid(mesh.nodes);
    elements = mesh.elements;
    edges = mesh.edges;
    domain = mesh.domain;
    
    % Save mesh data to .mat files
    save('mesh_nodes.mat', 'nodes');
    save('mesh_edges.mat', 'edges');
    save('mesh_elements.mat', 'elements');
    save('mesh_domain.mat', 'domain');
    
    % Call the Python function using pyrunfile
    result = pyrunfile("calc_saddle_point.py", "p_h_matrix, g_i_h_matrix");

    % Convert Python results to MATLAB arrays
    p_h_matrix = double(result.p_h_matrix);
    g_i_h_matrix = double(result.g_i_h_matrix);
end
