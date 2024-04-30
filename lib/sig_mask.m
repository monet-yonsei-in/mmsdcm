
function y_mask = sig_mask(obs_dt, dt, y_sig)
step = obs_dt/dt; % it should be integer 
t_mask = [1:step:size(y_sig,1)];
y_mask = y_sig(t_mask,:);
end
