function [pE,pC,M]=mnet_dcm_mmod_set_inital_values(mC,Aw,Cw,nobs,pF,obsmodel)
% [pE,pC,M]=mnet_dcm_mmod_set_inital_values(mC,Aw,Cw,nobs,pF,obsmodel)
% input: 
% mC    : model parameters
% nobs  : number of observation -> nbobs.vsdi, nobs.cai, nobs.bold
% pF    : fixed parameters
% Y     : observation signals
% U     : stimulus
% obsmodel : name of observation model

% output:
% pE    : prior mean 
% pC    : prior covariance
% M     : structures of parameters and conditions
%
% last updated Jiyoung Kang: 2019-09-13

if nargin < 1, mC = []; end
if nargin < 2, pF = []; end
if nargin < 3, Aw = []; end
if nargin < 4, Cw = []; end
if nargin < 5, nobs = []; end

% number of obs. signal
switch obsmodel
    case 'VSDI'
        %nobs2   =   nobs.vsdi;
        nobs2   =   sum(unique(mC.obs.vsdi)>0);
    case 'CaI'
        %nobs2   =   nobs.cai; 
        nobs2   =   sum(unique(mC.obs.cai)>0);
    case 'BOLD'
        %nobs2   =   nobs.bold; 
        nobs2   =   sum(unique(mC.obs.bold)>0);
    case 'VSDI_CaI'
        %nobs2   =   nobs.vsdi+nobs.cai; 
        nobs2   =   sum(unique(mC.obs.vsdi)>0)+sum(unique(mC.obs.cai)>0);
    case 'VSDI_BOLD'
        %nobs2   =   nobs.vsdi+nobs.bold; 
        nobs2   =   sum(unique(mC.obs.vsdi)>0)+sum(unique(mC.obs.bold)>0);
    case 'CaI_BOLD'
        %nobs2   =   nobs.cai+nobs.bold;  
        nobs2   =   sum(unique(mC.obs.cai)>0)+sum(unique(mC.obs.bold)>0);
    case 'VSDI_CaI_BOLD' 
        %nobs2   =   nobs.vsdi+nobs.cai+nobs.bold;    
        nobs2   =   sum(unique(mC.obs.vsdi)>0)+sum(unique(mC.obs.cai)>0)+sum(unique(mC.obs.bold)>0);
    otherwise
        nobs2=nobs;
end

[ns,~]  = size(Aw);     % number of neural populations

%% complete mC
mC.nobs = nobs; % it includes nobs.vsdi,nobs.cai,nobs.bold ...
% strength of stimulus
mC.C0         = 0.25; %20;

% vsdi
mC.Alpha      = 0.01 * ones(nobs.vsdi,1); % %[-] empirical setting % chk-H-related

% bold
mC.Hb_decay   = ones(nobs.bold,1) .* 0.64;
mC.Hb_transit = ones(nobs.bold,1) .* 2.0;

mC.Hb_epsilon = ones(1,1) .* 1.0;

mC.J_hd_exc   = ones(ns,1)  .* 0.1; % chk! %chg
mC.J_hd_inh   = ones(ns,1)  .* 0.1; % chk! %chg
mC.J_hd_ext   = ones(ns,1)  .* 0.1; % chk! %chg % nobs.blod?

%% complete pF
% If you want to free only for A,C use this default settting.
% You can also set pF in previouse stage of this function. 
if isempty(pF)
    % specify fixed parameters
    % ------------------------
    pF       = [];
    pF.R     = ones(ns,1) .* mC.R;  %[v]
    pF.H     = ones(ns,1) .* mC.H;  %[v]
    pF.T     = ones(ns,1) .* mC.T; 
    pF.Fr    = ones(ns,1) .* mC.Fr; %[v]
    pF.C0    = ones(1,1) * mC.C0;
    % for exponential A : A0*exp(A)
    % mC.A0   = 1; 
    pF.A0   = ones(1,1) * mC.A0;

    % VSDI
    pF.Alpha = ones(nobs.vsdi,1) .* mC.Alpha; %[v]
    pF.La    = ones(ns,1) .* mC.La; %[-]

    % CaI
    pF.Ca    = ones(1,1)   .* mC.Ca;      
    pF.KCa   = ones(1,1)   .* mC.KCa;      
    pF.TCa   = ones(1,1)   .* mC.TCa;      
    pF.RCa   = ones(ns,1)  .* mC.RCa; % Ca-dynamics is handdled at the neural population level  
    pF.GCa   = ones(1,1)   .* mC.GCa;      
    pF.ECa   = ones(1,1)   .* mC.ECa;      
    pF.thrCa = ones(1,1)   .* mC.thrCa;      

    pF.KF   = ones(1,1)    .* mC.KF;
    pF.Kd   = ones(1,1)    .* mC.Kd;
    pF.dF   = ones(1,1)    .* mC.dF;

    % BOLD
    pF.Hb_decay   = ones(1,1)   .* mC.Hb_decay;
    pF.Hb_transit = ones(1,1)   .* mC.Hb_transit;

    pF.Hb_epsilon = ones(1,1)   .* mC.Hb_epsilon;

    pF.J_hd_exc   = ones(ns,1)  .* mC.J_hd_exc; %ones(nobs.bold,1)  .* mC.J_hd_exc;
    pF.J_hd_inh   = ones(ns,1)  .* mC.J_hd_inh;
    pF.J_hd_ext   = ones(ns,1)  .* mC.J_hd_ext;

end

%%
% specify priors and initial states
% ---------------------------------
[pE,pC]      = mnet_dcm_mmod_priors(mC.A,mC.C,nobs,obsmodel);               % initial expectations and covariances of priors
[pE,pC]      = mnet_dcm_remove_prior(pE,pC,pF);    % remove fields of fixed-parameters from priors
[x,f,~,hx,fH]   = mnet_dcm_x_neural(pE,obsmodel,mC);     % initial states and equations of motion % x: initial states of neural pop. h: initial states of hidden components % fH: functions for dynamics of obs.

% intrinsic connectivity (ground truth; GT)
if ~isempty(Aw), pE.A = Aw; end
if ~isempty(Cw), pE.C = Cw; end

% specify M
% ---------
M         = [];
M.IS      = 'mnet_int_ode';
M.f       = f;
M.x       = x;

M.n       = length(spm_vec(x));           % the number of states
M.l       = ns;                           % the number of outputs + hidden states

M.nobs    = nobs2;                        % the number of outputs % nobs.vsdi+nobs.cai+nobs.bold
M.pE      = pE;
M.pC      = pC;
M.pF      = pF;
M.mC      = mC;
M.hE      = 6;                            % only for parameter estimation
M.hC      = 1/64;                         % only for parameter estimation
M.Nmax    = 1024;                         % only for parameter estimation
% M.delays  = ones(1,M.l) * 1;              % not used


switch obsmodel
    case{'VSDI'}        
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.vsdi  = 'mnet_dcm_gx_vsdi';                
        
        % fx functions for dynamics of obs. variables,
        M.fH.x = ''; % dummy to prevent error

    case{'CaI'}        
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.cai   = 'mnet_dcm_gx_cai';

        % hx: fx functions for dynamics of obs. variables,
        M.fH.h_cai   = fH.h_cai;
        
        % initial values for ca dyanmics and haemodynamics
        M.hx.x_cai  = hx.x_cai;

    case{'BOLD'}
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.bold  = 'mnet_dcm_gx_bold';

        % hx: fx functions for dynamics of obs. variables,
        M.fH.h_bold  = fH.h_bold;
        
        % initial values for ca dyanmics and haemodynamics
        M.hx.x_bold = hx.x_bold;

    case{'VSDI_CaI'}        
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.vsdi  = 'mnet_dcm_gx_vsdi';                
        M.hg.cai   = 'mnet_dcm_gx_cai';        

        % hx: fx functions for dynamics of obs. variables,
        M.fH.h_cai   = fH.h_cai;
        
        % initial values for ca dyanmics and haemodynamics
        M.hx.x_cai  = hx.x_cai;        

    case{'VSDI_BOLD'}        
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.vsdi  = 'mnet_dcm_gx_vsdi';                
        M.hg.bold  = 'mnet_dcm_gx_bold';

        % hx: fx functions for dynamics of obs. variables,
        M.fH.h_bold  = fH.h_bold;
        
        % initial values for ca dyanmics and haemodynamics
        M.hx.x_bold = hx.x_bold;

    case{'CaI_BOLD'}
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.cai   = 'mnet_dcm_gx_cai';
        M.hg.bold  = 'mnet_dcm_gx_bold';

        % hx: fx functions for dynamics of obs. variables,
        M.fH.h_cai   = fH.h_cai;
        M.fH.h_bold  = fH.h_bold;
        
        % initial values for ca dyanmics and haemodynamics
        M.hx.x_cai  = hx.x_cai;
        M.hx.x_bold = hx.x_bold;
        
        
    case{'VSDI_CaI_BOLD'}        
        M.IS       = 'mnet_int_ode_mmd';
        
        % gx for each obs. signals
        M.hg.vsdi  = 'mnet_dcm_gx_vsdi';                
        M.hg.cai   = 'mnet_dcm_gx_cai';
        M.hg.bold  = 'mnet_dcm_gx_bold';

        % hx: fx functions for dynamics of obs. variables,
        M.fH.h_cai   = fH.h_cai;
        M.fH.h_bold  = fH.h_bold;
        
        % initial values for ca dyanmics and haemodynamics
        M.hx.x_cai  = hx.x_cai;
        M.hx.x_bold = hx.x_bold;
end
end
