

%[cell_list, cell_list_by_idx, xlist,tlist]=get_node_list();
load("cell_list_from_y_0.1.mat");

N_cell = size(cell_list,1);

task_start = 1
NG_task_ids = [];
for k = task_start:N_cell

  fprintf("========================================== [ Cell No.: %d (Total %d) ] ========================================== \n", k, N_cell);

  r_cell = cell_list(k,:);
  try
      validate_region_cell(r_cell, k);
  catch ME
    fprintf("Cell No.: %s\n", k);
    fprintf("An error occurred: %s\n", ME.message);
    NG_task_ids = [NG_task_ids,k];

    save("NG_task_ids.mat","NG_task_ids");

  end

end
