
function get_best_sol(SimName)

f_mat = sprintf('%s.mat',SimName);
load(f_mat,'results','boption')

boption_new = extract_parm(results,boption);


end