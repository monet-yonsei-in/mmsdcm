function [pE,pC] = mnet_dcm_mmod_priors(A,C,nobs,model)
%med_vsdi_prior : prior for cai, vsdi, vsdi_cai
%   [pE,pC,x] = med_vsdi_priors(A,C,nobs)
%
%   inputs (binary constraint)
%       A           - intrinsic connectivity            (ns x ns), ns: the number of sources
%       C           - driving inputs                    (ns x nm), ni: the number of modulations
%       nobs        - number of observation
%       model       - name of model
%
%   outputs
%       pE (prior expectation) - f(x,u,P,M)
%           % connectivity parameters
%           % -----------------------
%           pE.A    - intrinsic connectivity                    (ns x ns)
%           pE.C    - driving inputs                            (ns x nm)
%
%           % synaptic parameters
%           % -------------------
%           pE.H    - maximum postsynaptic potential            (ns x 1)
%           pE.K    - inverse membrane time constanct for PSP   (ns x 1)
%           pE.R    - slope of sigmoid for activation function  (ns x 1)
%           pE.Fr   - maximal firing rate                       (ns x 1)
%
%       pC (prior covariances)
%
%   Jiyoung Kang, 2019-09-06

ns              = size(A,1); % the number of sources
ni              = size(C,2); % the number of external inputs (driving inputs)
%nobs            % the numebr of output

% prior expectations - f(x,u,P,M)
% -------------------------------
pE   = [];
if ~exist('A','var'),   pE.A = ones(ns,ns)*0;   else,   pE.A = double(A~=0)*0;   end
if ~exist('C','var'),   pE.C = ones(ns,ni)*0;   else,   pE.C = double(C~=0)*0;   end
pE.H  = ones(ns,1)  * 0;        % maximal postsynaptic potential (mV)
pE.T  = ones(ns,1)  * 0;        % time constanct for PSP (s)
pE.R  = ones(ns,1)  * 0;        % slope of sigmoid nueral activity
pE.Fr = ones(ns,1)  * 0;        % maximal firing rate (Hz)
pE.C0 = ones(1,1)   * 0;        % Coefficient of external input 

pE.A0 = ones(1,1)   * 0;

% prior covariances
% -----------------
pC    = [];
d1=32;
d2=64;
d3=128;
d4=256;

pC.A  = double(A~=0) / d1;      
pC.C  = double(C~=0) / d1;        
pC.H  = ones(ns,1)   / d4;      
pC.T  = ones(ns,1)   / d4;      
pC.R  = ones(ns,1)   / d4;      
pC.Fr = ones(ns,1)   / d4;      
pC.C0 = ones(1,1)    / d4;

% for A=A0*exp(A)
%   pE.A0 = ones(1,1)   * 0;        % Coefficient of effect size
pC.A0 = ones(1,1)    / d1;


if strcmp(model,'VSDI'), flag=[1 0 0]; end
if strcmp(model,'CaI'),  flag=[0 1 0]; end
if strcmp(model,'BOLD'), flag=[0 0 1]; end

if strcmp(model,'VSDI_CaI'),  flag=[1 1 0]; end
if strcmp(model,'VSDI_BOLD'), flag=[1 0 1]; end
if strcmp(model,'CaI_BOLD'),  flag=[1 0 1]; end

if strcmp(model,'VSDI_CaI_BOLD'),  flag=[1 1 1]; end

if(flag(1)==1)
% VSDI
   pE.La    = ones(ns,1)          * 0;     % lambda for VSDI
   pE.Alpha = ones(nobs.vsdi,1)   * 0;     % scaling factor for VSDI signal

   pC.La    = ones(ns,1)        / d1;
   pC.Alpha = ones(nobs.vsdi,1) / d1; 

end

if(flag(2)==1)
deno=d4;
% CaI
   % kinetics fx
   pE.Ca     = ones(1,1)   * 0;        % baseline calcium concentration
   pE.RCa    = ones(ns,1)  * 0;        % slope of sigmoid for Ca   
   pE.GCa    = ones(1,1)   * 0;        % maximal conductance for calcium ions   
   pE.KCa    = ones(1,1)   * 0;        % converting calcium current with scaling
   pE.TCa    = ones(1,1)   * 0;        % time constanct for calcium kernel
   pE.ECa    = ones(1,1)   * 0;        % calcium concentrations at rest
   pE.thrCa  = ones(1,1)   * 0;        % threshold for calcium dynamics

   % observation model gx   
   pE.KF = ones(1,1)   * 0;          % scale parameter for the fluorescence trace (i.e., maximal calcium signal)
   pE.Kd = ones(1,1)   * 0;          % dissociation constant (Yasuda et al., 2004, Sci STKE 2004: pl5.)
   pE.dF = ones(1,1)   * 0;          % offset parameter for the fluorescence trace
   
   pC.Ca     = ones(1,1)  / deno;      
   pC.RCa    = ones(ns,1) / deno;      
   pC.GCa    = ones(1,1)  / deno;      
   pC.KCa    = ones(1,1)  / deno;      
   pC.TCa    = ones(1,1)  / deno;      
   pC.ECa    = ones(1,1)  / deno;      
   pC.thrCa  = ones(1,1)  / deno;        % threshold for calcium dynamics
   
   pC.KF     = ones(1,1)  / deno;
   pC.Kd     = ones(1,1)  / deno;
   pC.dF     = ones(1,1)  / deno;

end

if(flag(3)==1) 
deno=d4;
% BOLD
   pE.Hb_decay   = ones(nobs.bold,1)   * 0;        % fx, signal decay  !
   pE.Hb_transit = ones(nobs.bold,1)   * 0;        % fx, transit time  !
   
   pE.J_hd_exc   = ones(ns,1)  * 0;        % coef. of excitation input on neural population i
   pE.J_hd_inh   = ones(ns,1)  * 0;        % coef. of inhibition input on neural population i
   pE.J_hd_ext   = ones(ns,1)  * 0;        % coef. of external   input on neural population i

   pE.Hb_epsilon = ones(1,1)   * 0;        % gx, region-specific ratios of intra- to extra-vascular signal !

   %
   pC.Hb_decay   = ones(nobs.bold,1)    / deno;
   pC.Hb_transit = ones(nobs.bold,1)    / deno;
   
   pC.J_hd_exc   = ones(ns,1)   / deno;
   pC.J_hd_inh   = ones(ns,1)   / deno;
   pC.J_hd_ext   = ones(ns,1)   / deno;
   
   pC.Hb_epsilon = ones(1,1)    / deno;

end


return;

end
