function [f,J,D] = mnet_dcm_hfx_bold(x,u,P,M)
% copied spm_fx_hdm(x,u,P,M)
%
% state equation for the hemodynamic model
% FORMAT [f] = spm_fx_hdm(x,u,P,M)
% x      - state vector
%   x(1) - vascular signal                                    s
%   x(2) - rCBF                                           log(f)
%   x(3) - venous volume                                  log(v)
%   x(4) - dHb                                            log(q)
% u      - input (neuronal activity)                      (u)
% P      - free parameter vector
%   P(1) - signal decay                                   d(ds/dt)/ds)
%   P(2) - autoregulation                                 d(ds/dt)/df)
%   P(3) - transit time                                   (t0)
%   P(4) - exponent for Fout(v)                           (alpha)
%   P(5) - resting oxygen extraction                      (E0)
%   P(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal   
%
%   P(6 + 1:m)   - input efficacies                       d(ds/dt)/du)
%
% y      - dx/dt
%__________________________________________________________________________
%
% Ref Buxton RB, Wong EC & Frank LR. Dynamics of blood flow and oxygenation
% changes during brain activation: The Balloon model. MRM 39:855-864 (1998)
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_fx_hdm.m 6856 2016-08-10 17:55:05Z karl $


% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%--------------------------------------------------------------------------

if isfield(P,'Hb_decay'),Hb_decay = exp(P.Hb_decay) * M.mC.Hb_decay;else,Hb_decay = M.pF.Hb_decay;end 
if isfield(P,'Hb_transit'),Hb_transit = exp(P.Hb_transit) * M.mC.Hb_transit;else,Hb_transit = M.pF.Hb_transit;end 

x       = spm_unvec(x,M.hx.x_bold); 


%if isstruct(P)
    Hb     = [0.64 0.32 2.00 0.32 0.4 1 1];
%     for i = 1:numel(P.decay)
%         H(1)   = H(1)*exp(P.decay(i));
%         H(3)   = H(3)*exp(P.transit(i));
%         f(i,:) = spm_fx_hdm(x(i,:),u(i),H);
%     end
%     f     = f(:);
%     return
%end


% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,2:end) = exp(x(:,2:end)); 

% signal decay
%--------------------------------------------------------------------------
sd       = Hb(1)*exp(Hb_decay);

% transit time
%--------------------------------------------------------------------------
tt       = Hb(3)*exp(Hb_transit);

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv       = x(:,3).^(1/Hb(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff       = (1 - (1 - Hb(5)).^(1./x(:,2)))/Hb(5);

% implement differential state equations
%--------------------------------------------------------------------------
%f(:,1)     = P(7:end)'*u(:) - P(1).* x(:,1) - P(2)*(x(:,2) - 1);
f(:,1)     = Hb(7:end)' * u' - sd .* x(:,1) - Hb(2)*(x(:,2) - 1); % n_g
f(:,2)     = x(:,1) ./ x(:,2);
f(:,3)     = (x(:,2) - fv) ./ (tt .* x(:,3));
f(:,4)     = (ff .* x(:,2) - fv .* x(:,4) ./ x(:,3)) ./(tt .* x(:,4));
%f        = f(:);

% % vascular signal
% f(:,4)   = s - sd.*x(:,4) - H(2)*(x(:,5) - 1); % n_g
% f(:,5)   = x(:,4)./x(:,5);                          %
% f(:,6)   = (x(:,5) - fv)./(tt.*x(:,6));
% f(:,7)   = (ff.*x(:,5) - fv.*x(:,7)./x(:,6))./(tt.*x(:,7));

f       = spm_vec(f);

% adjust motion for DEM (that uses time-bins as units of time)
%--------------------------------------------------------------------------
% try, global dt, f  = f*dt; end

return

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
