
function minDCM = opt1_3c1v(SimName0, run_options)

SimName = sprintf('%s_%s',SimName0, run_options.type);
f_log = sprintf('%s.txt',SimName);
f_mat = sprintf('%s.mat',SimName);
diary(f_log)
diary on 

GTfile = './GT2.mat';
nodedef = 'node_region1.csv';      % node
netdef = 'connect_region1.csv';    % connectivity
[M, Y, U, options] = prepare_dcm_input(GTfile, nodedef, netdef);

M.pE.A = zeros(size(M.pE.A));
M.pE.C = zeros(size(M.pE.C));
M.pE.T = zeros(size(M.pE.T));

calinfo.Y = Y;
calinfo.M = M;
calinfo.U = U;
calinfo.options = options;


parmchg = {'A(2,1)', 'A(3,1)', 'A(1,2)', 'A(3,2)', 'A(1,3)', 'A(2,3)', ... 
    'C(1)', 'C(3)', 'T(1)', 'T(2)', 'T(3)' };
for ii=1:length(parmchg)
    boption(ii).pEname= parmchg{ii}; 
    
    tmpname = strrep(parmchg{ii},')','');
    tmpname = strrep(tmpname,'(','');
    tmpname = strrep(tmpname,',','_');
    boption(ii).name =  tmpname; %['p' num2str(ii)];   

    if strncmp(parmchg{ii},'A',1)
        boption(ii).range = [-0.5 1.5 ]; 
    elseif strncmp(parmchg{ii},'C',1)
        boption(ii).range = [0.3 0.8 ];
    elseif strncmp(parmchg{ii},'T',1)
        boption(ii).range = [0.0 0.5 ];
    else
        disp(['chk-parmchg'])
    end
end

if strcmp(run_options.type,'fine')        
    % set BayesOpt range for finetuning
    f_mat1 = sprintf('%s_%s.mat',SimName0,'coarse');
    load(f_mat1,'results')    
    boption_new = extract_parm(results,boption);
    boption = boption_new;
end

results = mnet_dcm_search_param(calinfo, boption, run_options);
save(f_mat)

if strcmp(run_options.type,'fine')        
    % Get Best Results
    xvar=results.bestPoint;
    save_option = true;
    save_fname = sprintf('%s_DCM',SimName0);
    calinfo.boption = boption;
    calinfo.DCM_NMax = run_options.DCM_NMax;
    mF = mnet_dcm_cost(xvar, calinfo, save_option, save_fname);
    fprintf('%s F= %5.1f \n', SimName, -mF);
end

minDCM = -results.MinObjective;
fprintf('min_DCM %s F= %5.1f \n', SimName, minDCM); % for chk
diary off
end


function [results] = mnet_dcm_search_param(calinfo, boption, run_options)

for ii = 1:length(boption)
    
    name = boption(ii).name;
    range = boption(ii).range;
    cmd = sprintf("%s = optimizableVariable('%s',range);", name, name);
    eval(cmd);
    
    cmd = sprintf(" vars(%d) = %s ;", ii, name);
    eval(cmd);
end

calinfo.boption = boption;
calinfo.DCM_NMax = run_options.DCM_NMax;
fun =  @(x)mnet_dcm_cost(x, calinfo);
 
if run_options.do_parallel    
    results = bayesopt(fun, vars, 'UseParallel', true, 'MaxObjectiveEvaluations',run_options.Bayes_NMax);
else
    results = bayesopt(fun, vars, 'MaxObjectiveEvaluations',run_options.Bayes_NMax);
end

end

function mF = mnet_dcm_cost(xvar, calinfo, save_option, save_fname)
if nargin < 4, save_fname = []; end
if nargin < 3, save_option = false; end

bopt_test = [];
for i=1:length(xvar.Properties.VariableNames)
    eval(sprintf('bopt_test.%s=xvar.(xvar.Properties.VariableNames{%d});',xvar.Properties.VariableNames{i},i));
end

% hit
DCM_hit = estimate_DCM(bopt_test, calinfo);

mF = -1.0 * (DCM_hit.F);

if save_option 
    name = save_fname; %['DCM' num2str(-mF) ];
    save([name '.mat'])
    cmd=sprintf('print -dpng -r300 %s.png',name);
    eval(cmd);    
end    
% if( mF < -500.0 )
%     name = ['DCM' num2str(-mF) ];
%     save([name '.mat'])
%     cmd=sprintf('print -dpng -r300 %s.png',name);
%     eval(cmd);
% end
end



function [DCM] = estimate_DCM(bopt_test, calinfo)

M = calinfo.M;
Y = calinfo.Y;
U = calinfo.U;

M.pE.A = (M.pC.A ~= 0)*1.0;
M.pE.C = (M.pC.C ~= 0)*1.0;
M.pE.T = (M.pC.T ~= 0)*1.0;

for ii=1:size(calinfo.boption,2)
    pEname = calinfo.boption(ii).pEname;
    cmd = sprintf('M.pE.%s = bopt_test.%s;',pEname,calinfo.boption(ii).name);
    eval(cmd)
end

M.Nmax = calinfo.DCM_NMax;
[Ep, Cp, Eh, F, L] = mnet_spm_nlsi_GN(M, U, Y);

clear P
%P.name = name;
P.M = M;
P.Y = Y;
P.U = U;

P.Ep = Ep;                   % conditional expectation    E{P|y}
P.Cp = Cp;                   % conditional covariance     Cov{P|y}
P.Eh = Eh;                   % conditional log-precisions E{h|y}
P.Ce = exp(-Eh);             % ReML error covariance

% log evidence
% ------------
P.F = F;                    % Laplace log evidence
P.L = L;
P.ID = spm_data_id(P.Y.y);   % data ID
DCM = P;

end


function [M, Y, U, options] = prepare_dcm_input(GTfile, nodedef, netdef)

%% input signals
load(GTfile, 'GT2')
ycai = GT2.Y.ymask.cai(:, 1:3);
yvsdi = GT2.Y.ymask.vsdi(1:100:end, 1);
yhat = [yvsdi(:); ycai(:)];
Y.y = yhat(:);
Y.dt = 0.1;

U = GT2.U;
U.u(82:end) = [];
U.u(12:end) = 0;
% added for multistim
U.u = [U.u; U.u*0.8; U.u*1.2];
U.dt = 0.1;

% Read topology parameters & get dimension
[mC, Aw, Cw, nobs] = mnet_dcm_mmod_load_parameters(nodedef, netdef, U);

% specify obsmodel % updated to be able to apply for all cases
obsmodel = 'VSDI_CaI_BOLD';

pF = GT2.M.pF;
pF.R = pF.R(1:3);
pF.Fr = pF.Fr(1:3);
pF.H = pF.H(1:3);
pF.RCa = pF.RCa(1:3);

[pE, pC, M] = mnet_dcm_mmod_set_inital_values(mC, Aw, Cw, nobs, pF, obsmodel);
M.m = size(U.u,2);                   % the number of inputs
M.ns = size(U.u,1);                  % the number of samples
M.Hz = 1/U.dt;                       % sampling rate % delete Hz or dt
M.dt = U.dt;                         % delta t       %

M.obs = GT2.M.obs;
M.obs.mask = true; %false;  %if(isfield(M.obs,'mask')), M.obs = rmfield(M.obs,'mask'); end

M.obs.vsdi.dt = 0.1; % down sampled

options.obsmodel = obsmodel;
options.pE_load = 0;

options.name = 'bopt_cal';
options.used_load_prior = 1; % zero prior
options.flag_prt_fig = 1;
M.pE = GT2.Ep;

%M.hg = rmfield(M.hg,'vsdi');
M.hg = rmfield(M.hg,'bold');
%M.obs = rmfield(M.obs,'vsdi');
M.obs = rmfield(M.obs,'bold');
M.fH = rmfield(M.fH,'h_bold');

M.pF.KCa = M.mC.KCa * exp(M.pE.KCa);
M.pF.TCa = M.mC.TCa * exp(M.pE.TCa);

M.pE = rmfield(M.pE, 'KCa');
M.pE = rmfield(M.pE, 'TCa');
M.pC = rmfield(M.pC, 'KCa');
M.pC = rmfield(M.pC, 'TCa');


end


