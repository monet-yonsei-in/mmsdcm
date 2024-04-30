function [y] = mnet_dcm_gx_bold(x,u,P,M)
% Multimodal observation data generation: vsdi, cai, bold 
%   [y] = mnet_dcm_gx_vsdi_cai_bold(x,u,P,M)
%
%   x      - state vector
%   bold
%   x(:,1) - vascular signal                          s  
%   x(:,2) - rCBF                                  ln(f)
%   x(:,3) - venous volume                         ln(v)
%   x(:,4) - deoyxHb                               ln(q)
%
%   y(:,nobs.bold)  - BOLD signals
%                                               
%   See also mnet_dcm_fx_vsdi_cai_bold, med_dcm_mmod_priors
%   Jiyoung Kang, 2019-09-06
%   See also function spm_gx_fmri
%   [g,dgdx] = spm_gx_fmri(x,u,P,M)
%__________________________________________________________________________
%

%% Obsevation model
if isfield(P,'La') && isfield(M,'mC')
    if isfield(M.mC,'La')
        La = exp(P.La) .* M.mC.La;
    end
elseif isfield(M,'pF') && isfield(M.pF,'La')
    La = M.pF.La;
end % lambda for vsdi
if isfield(P,'Alpha') && isfield(M,'mC')
    if isfield(M.mC,'Alpha'),
        Alpha = exp(P.Alpha) .* M.mC.Alpha;
    end    
elseif isfield(M,'pF') && isfield(M.pF,'Alpha')
    Alpha = M.pF.Alpha;
end % scale parameter for vsdi
% if isfield(P,'dF'),dF = exp(P.dF) .* M.mC.dF;else,dF = M.pF.dF;end % offset parameter

%% BOLD
if isfield(P,'Hb_epsilon'),Hb_epsilon = exp(P.Hb_epsilon) .* M.mC.Hb_epsilon;else,Hb_epsilon = M.pF.Hb_epsilon;end % 

%y  = zeros(1,max(M.mC.Popu.bold));

if ~isfield(M,'mH')
    M.mH = false(M.l,1);
    if isfield(M,'nh')
        if M.nh > 0
            M.mH(1:M.nh) = true;
        end
    end
end

%% Bold
% Biophysical constants for 1.5T
%==========================================================================
 
% time to echo (TE) (default 0.04 sec)
%--------------------------------------------------------------------------
TE  = 0.04;
 
% resting venous volume (%)
%--------------------------------------------------------------------------
V0  = 4;

% estimated region-specific ratios of intra- to extra-vascular signal 
%--------------------------------------------------------------------------
ep  = exp(Hb_epsilon);
 
% slope r0 of intravascular relaxation rate R_iv as a function of oxygen 
% saturation S:  R_iv = r0*[(1 - S)-(1 - S0)] (Hz)
%--------------------------------------------------------------------------
r0  = 25;
 
% frequency offset at the outer surface of magnetized vessels (Hz)
%--------------------------------------------------------------------------
nu0 = 40.3; 
 
% resting oxygen extraction fraction
%--------------------------------------------------------------------------
E0  = 0.4;
 
%-Coefficients in BOLD signal model
%==========================================================================
k1  = 4.3*nu0*E0*TE;
k2  = ep*r0*E0*TE;
k3  = 1 - ep;
 
%-Output equation of BOLD signal model
%==========================================================================
v   = exp(x(:,3)); % exp(venous volume)
q   = exp(x(:,4)); % exp(deoyxHb)

y   = V0*(k1 - k1.*q + k2 - k2.*q./v + k3 - k3.*v);

return;
end

