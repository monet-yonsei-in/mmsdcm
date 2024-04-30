function [M, options, U] = get_config(config_file, U, cal_option)

if (nargin < 1)
    error("Configuration json file was not loaded")
end

% load configuration file
config = jsondecode(fileread(config_file));

if (nargout > 2 && (~config.flag_gen_external_input))
    error("U was not set")    
end
    
try config.flag_prt; catch, config.flag_prt = false; end
try config.flag_prt_fig; catch, config.flag_prt_fig = false; end

% % set time
% dt = config.time.dt;   
% t = ([1:config.time.total_step]) * dt + config.time.start;
% 
% % generate external input U
% if (config.flag_gen_external_input)     
%     fname = config.gen_external_input.func_name;
%     option_external_input = config.gen_external_input.option;
%     U = feval(fname, t, dt, option_external_input);         
% end

%% specify model
% -------------
try 
    obsmodel = cal_option.model;
catch
    obsmodel = config.obsmodel; %'CaI'; % 'VSDI_CaI_BOLD'; % 'VSDI','CaI','BOLD','VSDI_CaI_BOLD','VSDI_CaI', 'VSDI_BOLD'
end

%% generate observation signals
nodedef = config.model.node_file_name;        % node
netdef = config.model.connect_file_name;      % connectivity

% Read topology parameters & get dimension
[mC, Aw, Cw, nobs] = mnet_dcm_mmod_load_parameters(nodedef, netdef, U);

% specify obsmodel
%pF = config.param.pF;
pF = [];
[pE, pC, M] = mnet_dcm_mmod_set_inital_values(mC, Aw, Cw, nobs, pF, obsmodel);
M.m = size(U.u,2);                   % the number of inputs
M.ns = size(U.u,1);                  % the number of samples
M.Hz = 1/U.dt;                       % sampling rate % delete Hz or dt
M.dt = U.dt;                         % delta t       % 

options = config;
options.obsmodel = obsmodel;

end
