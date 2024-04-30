function [f,J,D] = mnet_dcm_hfx_cai(x,u,P,M)
% x  : ca concentration
% u  : membrane potential 
% M.x_cad : initial Ca concentration
%
% neural mass model for multimodal systems of vsdi and cai 
% CaI dynamics - 
%   [f,J,D] = mnet_dcm_hfx_cai(x,u,P,M)
%
%   x       = state vector
%   xv(:,1)  - Ca dynamics 
%
%   u(:,m)  = input time-serise of membrane potentail 
%
%   P.Ca    = baseline of concentration of calcium ions             (a.u.)
%
%   P.G2    = maximal conductance for calcium ions                  (mS/cm^2)
%
%   M       - model
%
%   f       - dx(t)/dt  = f(x(t))
%   J       - df(t)/dx(t)
%   D       - delay operator dx(t)/dt = f(x(t - d))
%                                     = D(d)*f(x(t))
%                                               
%   See also med_vsdi_cai_priors
%   Jiyoung Kang, 2019-04-16

% get dimensions and configure state variables
%---------------------------------------------------------------------------------------------------
x       = spm_unvec(x,M.hx.x_cai);     % Ca states

% Parameters
%===================================================================================================
% CaI  
if isfield(P,'Ca'),Ca = exp(P.Ca) * M.mC.Ca;else,Ca = M.pF.Ca;end % baseline concentration of calcium ions
if isfield(P,'RCa'),RCa = exp(P.RCa) * M.mC.RCa;else,RCa = M.pF.RCa;end % slope of sigmoid for NMDAR
if isfield(P,'GCa'),GCa = exp(P.GCa) * M.mC.GCa;else,GCa = M.pF.GCa;end % maximal conductance for calcium ions
if isfield(P,'KCa'),KCa = exp(P.KCa) * M.mC.KCa;else,KCa = M.pF.KCa;end % converting calcium current with scaling
if isfield(P,'TCa'),TCa = exp(P.TCa) * M.mC.TCa;else,TCa = M.pF.TCa;end 
if isfield(P,'ECa'),ECa = exp(P.ECa) * M.mC.ECa;else,ECa = M.pF.ECa;end 
if isfield(P,'thrCa'),thrCa = exp(P.thrCa) * M.mC.thrCa;else,thrCa = M.pF.thrCa;end 

% dynamics of calcium signal (Rahmati et al., 2016, PLOS Comp. Biol.)
%---------------------------------------------------------------------------------------------------
%    f(:,1)  =  -  KCa .* GCa * (u(:,1) - ECa) ./ (1 + exp(-(u(:,1) - thrCa) .* RCa)) ...
%               - (x(:,1) - Ca) ./ TCa;                                    % d/dt [Ca2+]

   f(:,1)  =  -  KCa .* GCa * (u' - ECa) ./ (1 + exp(-(u' - thrCa) .* RCa)) ...
              - (x(:,1) - Ca) ./ TCa;                                    % d/dt [Ca2+]

f       = spm_vec(f);

if nargout < 2; return, end

% Jacobian
%===================================================================================================
if isfield(M,'x'), x = spm_vec(M.x); else,  x = sparse(M.n,1); end
if isfield(M,'u'), u = spm_vec(M.u); else,  u = sparse(M.m,1); end
J  = spm_diff(M.f,x,u,P,M,1); % Jacobian (J)

if nargout < 3; return, end

% delays
%===================================================================================================
% Delay differential equations can be integrated efficiently (but
% approximately) by absorbing the delay operator into the Jacobian
%
%    dx(t)/dt     = f(x(t - d))
%                 = Q(d)f(x(t))
%
%    J(d)         = Q(d)df/dx
%---------------------------------------------------------------------------------------------------
% Implement: dx(t)/dt = f(x(t - d)) = inv(1 + D.*dfdx)*f(x(t))
%                     = Q*f = Q*J*x(t)
%---------------------------------------------------------------------------------------------------
D  = spm_dcm_delay(P,M,J);
return;
