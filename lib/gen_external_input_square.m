function U = gen_external_input_square(t, dt, option)
% generate squre input

tstart = option.stim_start_nstep;
tend = option.stim_end_nstep;

u = zeros(size(t))';
u(tstart:tend) = 1;
U.u = u;
U.dt = dt;

end
