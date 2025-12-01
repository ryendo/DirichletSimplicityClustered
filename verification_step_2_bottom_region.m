
my_intlab_config_server

[cell_list, cell_list_by_idx, xlist,tlist]=get_node_list();

N_cell = size(cell_list,1);

for k = 1:N_cell
  disp("Task No.")
  disp([k,N_cell])

  r_cell = cell_list(k,:);
  validate_region_cell(r_cell, k);

end
