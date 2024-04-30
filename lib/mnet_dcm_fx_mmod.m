function [f,J,D] = mnet_dcm_fx_mmod(x,u,P,M,OPT)
% neural mass model for multimodal systems of vsdi and cai 
%
%   x       = state vector
%   x(:,1)  - cross-membrane potential difference               (mV)
%   x(:,2)  - cross-membrane current                            (mA)
%
% deleted
% %   x(:,3)  - concentration of calcium ions                     (a.u.)
% %   haemodynamics
% %   x(:,4)  - vascular signal                          s  
% %   x(:,5)  - rCBF                                  ln(f)
% %   x(:,6)  - venous volume                         ln(v)
% %   x(:,7)  - deoyxHb                               ln(q)
% end of deletion
%
%
%   u(:,m)  = extrinsic inputs 
%
%   P.H     = maximum value of response kernel 
%               - H1: maximal postsynaptic potential                (mV)
%
%   P.T     = time constants 
%
%   P.R     = slope of sigmoid activation function                  (H x gamma: connectivity A)
%
%   P.Fr    = maximal firing rate                                   (Hz)
%
%
%   M       - model
%
%   f       - dx(t)/dt  = f(x(t))
%   J       - df(t)/dx(t)
%   D       - delay operator dx(t)/dt = f(x(t - d))
%                                     = D(d)*f(x(t))
%
% When argv>4 
% FORMAT [u,v,w] = spm_fx_cmc_tfm(x,u,P,M,'activity')
% u  - intrinsic presynaptic input (inhibitory)
% v  - intrinsic presynaptic input (excitatory)
% w  - extrinsic presynaptic input
%
%                                               
%   See also spm_fx_cmc_tfm.m
%   
%   last updated 2019-09-11, Jiyoung Kang

% get dimensions and configure state variables
%---------------------------------------------------------------------------------------------------
x       = spm_unvec(x,M.x);                                       % the neuronal states

% Parameters
%===================================================================================================
if isfield(P,'H'), H = exp(P.H)   .* M.mC.H;else, H = M.pF.H;  end % maximal postsynaptic potential (mV)
if isfield(P,'T'), T = exp(P.T)   .* M.mC.T;else, T = M.pF.T;  end % inverse membrane time constanct for PSP
if isfield(P,'R'), R = exp(P.R)   .* M.mC.R;else, R = M.pF.R;  end % slope of sigmoid for AMPAR
if isfield(P,'Fr'),Fr = exp(P.Fr) .* M.mC.Fr;else,Fr = M.pF.Fr;end % maximal firing rate (Hz)
if isfield(P,'A0'),A0 = exp(P.A0) .* M.mC.A0;else,A0 = M.pF.A0;end % baseline of A
if isfield(P,'C0'),C0 = exp(P.C0) .* M.mC.C0;else,C0 = M.pF.C0;end % baseline of C

% intrinsic connectivity ReLU-A (Rectified Linear Unit A) or absolute A
% ---------------------------------------------------------------------
%A = A0 * abs(P.A);
%A0 = 0.0;
A  = exp(P.A) * A0;
A(M.pC.A == 0)  = 0;

%A(M.mH,:) = exp(P.A(M.mH,:)) * A0;
%A(M.pC.A == 0)  = 0;

% external input
% --------------
C               = exp(P.C) * C0;
C(M.pC.C == 0)  = 0;



% pre-synaptic inputs: s(V) x max. firing rate
%---------------------------------------------------------------------------------------------------
thr = M.mC.thr;
S       = 1 ./ (1 + exp(-(x(:,1) - thr) .* R)) .* Fr; %1 ./ (1 + exp(-(x(:,1) - M.mC.thr) .* R)) .* Fr;         % faster than sigmf(x,[A,C])
%S = 1./(1 + exp(-R * (x(:,1) - thr))) - 1./(1 + exp(R*thr)); % changed like erp mod. 2020-08-24 jykang
A = A .* repmat(M.mC.ExIn',M.l,1);

% Motion of states: f(x)
%===================================================================================================

% post-synaptic depolarization
%---------------------------------------------------------------------------------------------------
f(:,1)  = x(:,2);                                                   % d/dt v
depolV  = H .* (A * S);                                             % sum of depolarizations from neurons
depolV  = depolV + C .* u .* Fr .* H ; %C .* u;   %chk! mod. jykang:  %C .* u;  org: C .* u .* Fr .* H;  % depolarization
%f(:,2)  = (depolV - 2 * x(:,2) - (x(:,1) - M.x(:,1)) .* K) .* K;    % d/dt i
f(:,2)  = (depolV - 2 * x(:,2) - (x(:,1) - M.x(:,1)) ./T ) ./ T;    % d/dt i

f       = spm_vec(f);

if nargin > 4
    
    A_exc = A; A_exc(A_exc < 0) = 0;
    A_inh = A; A_inh(A_inh > 0) = 0;
%     u = J_hd_exc .* (A_exc * S); % intrinsic excitation
%     v = J_hd_inh .* (A_inh * S); % intrinsic inhibition
%     w = J_hd_ext .* C .* u;      % extrinsic excitation

    u2 = (A_exc * S); % intrinsic excitation
    v  = (A_inh * S); % intrinsic inhibition
    w  = C .* u;      % extrinsic excitation

    %s  = u + v + w    % vasular sig
   
    f = u2; J = v; D = w;         
    return
end


% % dynamics of calcium signal (Rahmati et al., 2016, PLOS Comp. Biol.)
% %---------------------------------------------------------------------------------------------------
% %    f(:,3)  =  -  KCa .* GCa * (x(:,1) - ECa) ./ (1 + exp(-(x(:,1) - M.thrCa(1)) .* RCa)) ...
% %               - (x(:,3) - Ca) ./ TCa;                                    % d/dt [Ca2+]
%  
% %           
% % implement differential state equation f = dx/dt (hemodynamic)
% % cf. spm_fx_fmri
% % x      - state vector
% %   x(:,1)Not used - excitatory neuronal activity            ue
% %   x(:,2)4 - vascular signal   h1                       s  
% %   x(:,3)5 - rCBF              h2                    ln(f)
% %   x(:,4)6 - venous volume     h3                    ln(v)
% %   x(:,5)7 - deoyxHb           h4                    ln(q)
% %  [x(:,6) - inhibitory neuronal activity             ui
% %                               h1-h4 notations : H.Wei et al. 2019
% %
% % hemodynamic parameters
% %--------------------------------------------------------------------------
% %   H(1) - signal decay                                   d(ds/dt)/ds)
% %   H(2) - autoregulation                                 d(ds/dt)/df)
% %   H(3) - transit time                                   (t0)
% %   H(4) - exponent for Fout(v)                           (alpha)
% %   H(5) - resting oxygen extraction                      (E0)
% %   H(6) - ratio of intra- to extra-vascular components   (epsilon)
% %          of the gradient echo signal
% %--------------------------------------------------------------------------
 
% A_exc = A; A_exc(A_exc < 0) = 0;
% A_inh = A; A_inh(A_inh > 0) = 0;
% s = J_hd_exc .* (A_exc * S) + J_hd_inh .* (A_inh * S);
% s  = s + J_hd_ext .* C .* u  ;   % external
% 
% % vascular signal
% f(:,4)   = s - sd.*x(:,4) - H(2)*(x(:,5) - 1); % n_g
% f(:,5)   = x(:,4)./x(:,5);                          %
% f(:,6)   = (x(:,5) - fv)./(tt.*x(:,6));
% f(:,7)   = (ff.*x(:,5) - fv.*x(:,7)./x(:,6))./(tt.*x(:,7));


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
