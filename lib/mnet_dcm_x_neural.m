function [x,f,h,hx,fH] = mnet_dcm_x_neural(P,model,mC)
% Returns the state and equation of neural mass models
% FORMAT [x,f,h] = spm_dcm_x_neural(P,'model',mC)
%
%  P      - parameter structure
% 'model' - 'ERP','SEP','CMC','LFP','CMM','NNM', 'MFM', 'CMM NMDA' or 'CASD'
%  mC     - constraint parameters in the model
%
% x   - initial states
% f   - state equation dxdt = f(x,u,P,M)  - synaptic activity
% h   - state equation dPdt = f(x,u,P,M)  - synaptic plasticity
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging
%
%   Kyesam Jung added 'CASD' model (27 June 2018).

% Karl Friston
% $Id: spm_dcm_x_neural.m 6855 2016-08-06 10:06:35Z karl $

% paramteric state equation
%--------------------------------------------------------------------------
h  = [];

% assemble initial states for generic models
%==========================================================================
if isstruct(model)
    
    P.A{1} = 1;
    for i = 1:numel(model)
        x{i} = ksj_dcm_x_neural(P,model(i).source);
    end
    
    % general (multi-model) equations of motion
    %----------------------------------------------------------------------
    f  = 'spm_fx_gen'; 
    return
    
end



% initial state and equation
%--------------------------------------------------------------------------
switch lower(model)
    
    % linear David et al model (linear in states)
    %======================================================================
    case{'erp'}
        
        % inital states and equations of motion
        %------------------------------------------------------------------
        n  = length(P.A{1});                          % number of sources
        m  = 9;                                       % number of states
        x  = sparse(n,m);
        
        f  = 'spm_fx_erp';
        
        
    % linear David et al model (linear in states) - fast version for SEPs
    %======================================================================
    case{'sep'}
        
        % inital states
        %------------------------------------------------------------------
        n  = length(P.A{1});                          % number of sources
        m  = 9;                                       % number of states
        x  = sparse(n,m);
        
        f  = 'spm_fx_sep';
        
    % Linear in states ? canonical microcircuit
    %======================================================================
    case{'cmc'}
        
        % inital states
        %------------------------------------------------------------------
        n  = length(P.A{1});                          % number of sources
        m  = 8;                                       % number of states
        x  = sparse(n,m);
        
        f  = 'spm_fx_cmc';
        %%% f  = 'spm_fx_cmc_2014'; %%%
        
    % Linear in states ? canonical microcircuit with plasticity
    %======================================================================
    case{'tfm'}
        
        % inital states
        %------------------------------------------------------------------
        n  = length(P.A{1});                          % number of sources
        m  = 8;                                       % number of states
        x  = sparse(n,m);
        
        f  = 'spm_fx_cmc_tfm';
        h  = 'spm_fp_cmc_tfm';
        
    % linear David et al model (linear in states) - with self-inhibition
    %======================================================================
    case{'lfp'}
        
        % inital states
        %------------------------------------------------------------------
        n  = length(P.A{1});                          % number of sources
        m  = 13;                                      % number of states
        x  = sparse(n,m);
        
        f  = 'spm_fx_lfp';
        
        
    % Neural mass model (nonlinear in states)
    %======================================================================
    case{'nmm'}
        
   
        % get initialisation from full mean-field model
        %------------------------------------------------------------------
        x  = spm_x_mfm(P);
        
        % remove dispersion and fix the covariance of the states (Cx)
        %--------------------------------------------------------------------------
        x  = x{1};
        f  = 'spm_fx_mfm';
        
            % Neural mass model (nonlinear in states)
    %======================================================================
    case{'nmda'}
        
   
        % get initialisation from full mean-field model
        %------------------------------------------------------------------
        x  = spm_x_nmda(P);
        
        % remove dispersion and fix the covariance of the states (Cx)
        %--------------------------------------------------------------------------
        x  = x{1};
        f  = 'spm_fx_nmda';
        
    % Canonical mass model (nonlinear in states)
    %======================================================================
    case{'cmm'}
        
        % inital states and model
        %------------------------------------------------------------------
        x  = spm_x_cmm(P);
        f  = 'spm_fx_cmm';
        
  % Canonical mass model with NMDA (nonlinear in states)
    %======================================================================
    case{'cmm_nmda'}
        
        % inital states and model
        %------------------------------------------------------------------
        x  = spm_x_cmm_NMDA(P);
        f  = 'spm_fx_cmm_NMDA';
        
        
    % Mean field model (nonlinear in states) - with covariance
    %======================================================================
    case{'mfm'}
        
        % inital states and model
        %------------------------------------------------------------------
        x  = spm_x_mfm(P);
        f  = 'spm_fx_mfm';
        
            % Mean field model (nonlinear in states) - with covariance
    %======================================================================
    case{'null'}
        
        % inital states and model
        %------------------------------------------------------------------
        x  = sparse(size(P.A,1),1);
        f  = 'spm_fx_null';
        
    case{'cai'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)
        
        hx.x_cai  = ones(mC.nobs.cai,1) * mC.Ca;  % % baseline concentration of calcium ions    (a.u.)

        f           = 'mnet_dcm_fx_mmod';
        fH.h_cai    = 'mnet_dcm_hfx_cai';
        
    case{'vsdi'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)

        hx.x_vsdi = [];     % no h_vsdi dynamics for vsdi

        f           = 'mnet_dcm_fx_mmod';
        fH.h_vsdi   = ''; % to prevent error

    case{'bold'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)
        
        hx.x_bold = zeros(mC.nobs.bold,4); % initial values for HbA-dynamics: vascular signal, rCBF, venous volume, deoyxHb

        f           = 'mnet_dcm_fx_mmod';
        fH.h_bold   = 'mnet_dcm_hfx_bold';
        
    case{'vsdi_cai'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)

        %hx.x_vsdi = [];     % no h_vsdi dynamics for vsdi
        hx.x_cai  = ones(mC.nobs.cai,1) * mC.Ca;  % % baseline concentration of calcium ions    (a.u.)

        f           = 'mnet_dcm_fx_mmod';
        %fH.h_vsdi   = ''; % to prevent error
        fH.h_cai    = 'mnet_dcm_hfx_cai';

    case{'vsdi_bold'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)
        
        hx.x_bold = zeros(mC.nobs.bold,4); % initial values for HbA-dynamics: vascular signal, rCBF, venous volume, deoyxHb

        f           = 'mnet_dcm_fx_mmod';
        fH.h_bold   = 'mnet_dcm_hfx_bold';
        
    case{'cai_bold'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)
        
        hx.x_cai  = ones(mC.nobs.cai,1) * mC.Ca;  % % baseline concentration of calcium ions    (a.u.)
        hx.x_bold = zeros(mC.nobs.bold,4); % initial values for HbA-dynamics: vascular signal, rCBF, venous volume, deoyxHb

        f           = 'mnet_dcm_fx_mmod';
        fH.h_cai    = 'mnet_dcm_hfx_cai';
        fH.h_bold   = 'mnet_dcm_hfx_bold';
        
    case{'vsdi_cai_bold'}
        x = zeros(size(P.A,1),2);
        x(:,1) = mC.Vm;     % cross-membrane potential difference       (mV)
        x(:,2) = 0;         % cross-membrane current                    (mA)
        
        hx.x_cai  = ones(mC.nobs.cai,1) * mC.Ca;  % % baseline concentration of calcium ions    (a.u.)
        hx.x_bold = zeros(mC.nobs.bold,4); % initial values for HbA-dynamics: vascular signal, rCBF, venous volume, deoyxHb

        f           = 'mnet_dcm_fx_mmod';
        fH.h_cai    = 'mnet_dcm_hfx_cai';
        fH.h_bold   = 'mnet_dcm_hfx_bold';
        
    otherwise
        warndlg('Unknown model')
end
