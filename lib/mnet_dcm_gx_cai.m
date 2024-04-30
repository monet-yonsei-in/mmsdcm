function [y] = mnet_dcm_gx_cai(x,u,P,M)
% Multimodal observation data generation: vsdi, cai, bold 
%   [y] = mnet_dcm_gx_vsdi_cai_bold(x,u,P,M)
%
%   x      - state vector
%   x(:,1) - CaI                           (mM)

%   y(:,nobs.CaI)         - CaI signals

%   See also mnet_dcm_fx_vsdi_cai_bold, med_dcm_mmod_priors
%   Jiyoung Kang, 2019-09-06
%   See also function spm_gx_fmri
%   [g,dgdx] = spm_gx_fmri(x,u,P,M)
%__________________________________________________________________________
%

%% Calcium sensitive dye signal of time (Rahmati et al., 2016, PLOS Comp. Biol.)
if isfield(P,'KF'),KF = exp(P.KF) * M.mC.KF;else,KF = M.pF.KF;end % scale parameter for the fluorescence trace
if isfield(P,'Kd'),Kd = exp(P.Kd) * M.mC.Kd;else,Kd = M.pF.Kd;end % dissociation constant (Yasuda, Svoboda, et al., 2004, Sci STKE 2004: pl5.)
if isfield(P,'dF'),dF = exp(P.dF) * M.mC.dF;else,dF = M.pF.dF;end % offset parameter for the fluorescence trace
% if isfield(P,'Ca'),Ca = exp(P.Ca) * M.mC.Ca;else,Ca = M.pF.Ca;end % baseline concentration of calcium ions
% if isfield(P,'R1'),R1 = exp(P.R1) * M.mC.R1;else,R1 = M.pF.R1;end % slope of sigmoid for AMPAR
% if isfield(P,'Fr'),Fr = exp(P.Fr) * M.mC.Fr;else,Fr = M.pF.Fr;end % maximal firing rate (Hz)

%% CaI
% Observation model (Rahmati et al., 2016, PLOS Comp. Biol.)
% ----------------------------------------------------------
y  = KF .* x(:,1) ./ (x(:,1) + Kd) + dF;

% Linearity (Yaksi and Friedrich, 2006, Nature Methods) - deconvolution
% ---------------------------------------------------------------------
% y  = 1 ./ (1 + exp(-(x(:,1) - M.thr(1)) .* R1)) .* Fr .* KF;

return;
end

