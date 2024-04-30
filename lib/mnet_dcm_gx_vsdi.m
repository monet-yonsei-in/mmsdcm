function [y] = mnet_dcm_gx_vsdi(x,u,P,M)
% Multimodal observation data generation: vsdi, cai, bold 
%   [y] = mnet_dcm_gx_vsdi_cai_bold(x,u,P,M)
%
%   x      - state vector
%   x(:,1) - voltage                       (mV)
%   x(:,2) - current                       (mA)
%
%   y(:,nobs.vsdi)                   - VSDI signals
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

% Observation model for VSDI
% --------------------------
y  = zeros(1,max(M.mC.Popu.vsdi));
dFvsd = zeros(1,max(M.mC.Popu.vsdi));
for nY = 1:max(M.mC.Popu.vsdi)
    lgc      = M.mC.Popu.vsdi == nY;
    %org: y(1,nY)  = sum(x(lgc,1) .* La(lgc));
    %added jykang to change La
    scale  = sum(La(lgc));
    y(1,nY)  = sum(x(lgc,1) .* La(lgc) /scale);
    dFvsd(1,nY) = sum(M.x(lgc,1) .* La(lgc)/scale); %sum(M.x(lgc,1) .* La(lgc) * R2);
end
% dF = mean(dF);
% y = 1 ./ (1 + exp(-(y - 0) .* R2)) .* M.yRng(2) + M.yRng(1); % faster than sigmf(x,[A,C])
y = Alpha' .* (y - dFvsd); %R2 .* y - dF;

% added to fit scale
if(isfield(M,'scale_vsdi2cai'))
    y=y*M.scale_vsdi2cai;
end

%%
%y=[y];


return;
end

