function [y,xstat,xh,y_bf_mask] = mnet_int_ode_mmd(P,M,U)
% integrates a MIMO nonlinear system (using classical ODE solvers) for
% multimodal dynamics of electrophysiological dataset
%
% [y,xstat,xh] = mnet_int_ode_mmd(P,M,U)
% integrates a MIMO nonlinear system (using classical ODE solvers)
%
% modified to calculate
%    dx/dt      = f(x,u,P,M)
%    dx_cai/dt  = h_cai(x_cai,u,P,M)   % for Ca-dynamics
%    dx_bold/dt = h_bold(x_bold,u,P,M)  % for haemodynamics
%    y          = g(x,u,P,M) or y     = g_vsdi(x,u,P,M)
%    y_cai      = g_cai(x_cai,u,P,M)
%    y_bold     = g_bold(x_bold,u,P,M) % u = vascular-signals
%
% input 
%    P   - model parameters
%    M   - model structure
%    U   - input structure or matrix
% 
% output
%    y      = g(x,u,P,M)  % multiple observation signals
%    y      = [y; y_cai; y_bold]
%    xstat(ns,1:2)            : neural states ns=number of neural populations 1=PSP, 2=current 
%    xh.x_cai(ns,1)           : Ca cuurent    ns=number of neural populations 
%    xh.x_bold(nobs.bold,1:4) : haemodynamics nobs.bold = number of regions for BOLD signals
%%
%__________________________________________________________________________
% Integrates the MIMO system described by
%
% using an Runge-Kutta(4,5) scheme over the times implicit in the input.
% ode45 is based on an explicit Runge-Kutta (4,5) formula, the Dormand-
% Prince pair. It is a one-step solver - in computing y(tn), it needs only 
% the solution at the immediately preceding time point y(tn-1). In general,
% ode45 is the best function to apply as a "first try" for most problems.
%
% ode113 is a variable order Adams-Bashforth-Moulton PECE solver. It may be
% more efficient than ode45 at stringent tolerances and when the ODE file 
% function is particularly expensive to evaluate. ode113 is a multi-step 
% solver - it normally needs the solutions at several preceding time points
% to compute the current solution
%
% see also ode45; ode113
%--------------------------------------------------------------------------
%
% SPM solvers or integrators
%
% spm_int_ode:  uses ode45 (or ode113) which are one and multi-step solvers
% respectively.  They can be used for any ODEs, where the Jacobian is
% unknown or difficult to compute; however, they may be slow.
%
% spm_int_J: uses an explicit Jacobian-based update scheme that preserves
% nonlinearities in the ODE: dx = (expm(dt*J) - I)*inv(J)*f.  If the
% equations of motion return J = df/dx, it will be used; otherwise it is
% evaluated numerically, using spm_diff at each time point.  This scheme is
% infallible but potentially slow, if the Jacobian is not available (calls
% spm_dx).
%
% spm_int_E: As for spm_int_J but uses the eigensystem of J(x(0)) to eschew
% matrix exponentials and inversion during the integration. It is probably
% the best compromise, if the Jacobian is not available explicitly.
%
% spm_int_B: As for spm_int_J but uses a first-order approximation to J
% based on J(x(t)) = J(x(0)) + dJdx*x(t).
%
% spm_int_U: like spm_int_J but only evaluates J when the input changes.
% This can be useful if input changes are sparse (e.g., boxcar functions).
% spm_int_U also has the facility to integrate delay differential equations
% if a delay operator is returned [f J D] = f(x,u,P,M)
%
% spm_int:   Fast integrator that uses a bilinear approximation to the 
% Jacobian evaluated using spm_bireduce. This routine will also allow for
% sparse sampling of the solution and delays in observing outputs
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_int_ode.m 5219 2013-01-29 17:07:07Z spm $

% last updated 2019-09-18 Jiyoung Kang

% 
xh = [];
y_bf_mask = [];

% convert U to U.u if necessary
%--------------------------------------------------------------------------
if ~isstruct(U),  U.u  = U; end
try, U.dt; catch, U.dt = 1; end

% sample times
%--------------------------------------------------------------------------
ns      = size(U.u,1);
tspan   = (1:ns)*U.dt;
 
% state equation; add [0] states if not specified
%--------------------------------------------------------------------------
try
    f   = fcnchk(M.f,'x','u','P','M');
catch
    f   = inline('sparse(0,1)','x','u','P','M');
    M.n = 0;
    M.x = sparse(0,0);
end

% added by jykang
if isfield(M.fH,'h_cai')  %if ~isfield(M.fH,'h_cai'), M.fH.h_cai=[]; end 
    % and CaI dynamics     
    %--------------------------------------------------------------------------
    try
        h_cai   = fcnchk(M.fH.h_cai,'x','u','P','M');
        g_cai   = fcnchk(M.hg.cai,'x','u','P','M');                    
    catch
        h_cai   = [];
        g_cai   = [];
    end
end
if isfield(M.fH,'h_bold')
    % and Haemodynamics 
    %--------------------------------------------------------------------------
    try
        h_bold   = fcnchk(M.fH.h_bold,'x','u','P','M');
        g_bold   = fcnchk(M.hg.bold,'x','u','P','M');                    
    catch
        h_bold   = [];
        g_bold   = [];
    end
end
% end of addtion, jykang

% and output nonlinearity
%--------------------------------------------------------------------------
try
    g   = fcnchk(M.g,'x','u','P','M');
catch
    try
        g   = fcnchk(M.hg.vsdi,'x','u','P','M');
    catch
        g   = [];
    end
end
 
%%
% Initial states and inputs
%--------------------------------------------------------------------------
try
    x   = M.x;
catch
    x   = sparse(0,1);
    M.x = x;
end

% added by jykang
if isfield(M.fH,'h_cai')  %if ~isfield(M.fH,'h_cai'), M.fH.h_cai=[]; end 
    % and CaI dynamics     
    %--------------------------------------------------------------------------
    try
        x_cai   = M.hx.x_cai;
    catch
        x_cai   = sparse(0,1);
        M.hx.x_cai = x_cai;
    end
end
if isfield(M.fH,'h_bold')
    % and Haemodynamics 
    %--------------------------------------------------------------------------
    try
        x_bold   = M.hx.x_bold;
    catch
        x_bold   = sparse(0,1);
        M.hx.x_bold = x_bold;
    end
end
% end of addition jykang


% ODE45 functional form (note, the third argument of the ODE function is
% used by ode15i (i.e, OPTIONS).
%--------------------------------------------------------------------------
ode   = inline('spm_vec(f(spm_unvec(x,M.x),U.u(ceil(t/U.dt),:),P,M))',...
               't','x','OPTIONS','P','M','U','f');
           

OPTIONS = odeset;
[t,x]   = ode113(ode,tspan,spm_vec(M.x),OPTIONS,P,M,U,f);

% Give the integrator something to work with for 'flat' inputs:
%--------------------------------------------------------------------------
if norm(x,1) < exp(-16)
   U.u   = spm_conv(U.u,8,0); 
   [t,x] = ode113(ode,tspan,spm_vec(M.x),OPTIONS,P,M,U,f);
end

% added by jykang
if isfield(M.fH,'h_cai')  
    % solve ODE for Ca-dynamics
    U_cai.u=x(:,1:M.l); % M.l = # of neural pop.
    U_cai.dt=U.dt;
    %--------------------------------------------------------------------------
    ode_cai   = inline('spm_vec(h_cai(spm_unvec(x,M.hx.x_cai),U_cai.u(ceil(t/U_cai.dt),:),P,M))',...
                   't','x','OPTIONS','P','M','U_cai','h_cai');

    OPTIONS = odeset;
    [t,x_cai]   = ode113(ode_cai,tspan,spm_vec(M.hx.x_cai),OPTIONS,P,M,U_cai,h_cai);
    
    % Give the integrator something to work with for 'flat' inputs:
    %--------------------------------------------------------------------------
    if norm(x_cai,1) < exp(-16)
       U.u   = spm_conv(U.u,8,0); 
       [t,x_cai]   = ode113(ode_cai,tspan,spm_vec(M.hx.x_cai),OPTIONS,P,M,U,h_cai);
    end    
end

if isfield(M.fH,'h_bold')  
    
    %y(:,i) = spm_vec(g(spm_unvec(x(i,:),M.x),U.u(i,:),P,M));
    t_max_step     = size(U.u,1);

    % other alghrithm could be better! - chk
    for t_tmp=1:t_max_step
        x_tmp       = x(t_tmp,:);
        u_tmp       = U.u(t_tmp,:);
        [p,q,r]     = f(x_tmp,u_tmp,P,M,'activity');
        s(:,t_tmp)  = M.mC.J_hd_exc .* p + M.mC.J_hd_inh .* q + M.mC.J_hd_ext .* r;           
%         tmp_p(:,t_tmp)=p';
%         tmp_q(:,t_tmp)=q';
%         tmp_r(:,t_tmp)=r';
    end        
    
    
    % vascular signal for each region
    % --------------------------
    s_region       = zeros(M.mC.nobs.bold,t_max_step);
    for n_Sregion = 1:(M.mC.nobs.bold)
        lgc      = (M.mC.Popu.bold == n_Sregion);
        for t_tmp=1:t_max_step
            s_region(n_Sregion,t_tmp)  = sum(s(lgc,t_tmp)); 
        end        
    end
    
    U_bold.u    = s_region';
    U_bold.dt   = U.dt;
    % solve ODE for BOLD-dynamics
    %--------------------------------------------------------------------------
    ode_bold   = inline('spm_vec(h_bold(spm_unvec(x,M.hx.x_bold),U_bold.u(ceil(t/U_bold.dt),:),P,M))',...
                   't','x','OPTIONS','P','M','U_bold','h_bold');

    OPTIONS = odeset;
    [t,x_bold]   = ode113(ode_bold,tspan,spm_vec(M.hx.x_bold),OPTIONS,P,M,U_bold,h_bold);
    
    % Give the integrator something to work with for 'flat' inputs:
    %--------------------------------------------------------------------------
    if norm(x_bold,1) < exp(-16)
       U.u   = spm_conv(U.u,8,0); 
       [t,x_bold]   = ode113(ode_bold,tspan,spm_vec(M.hx.x_bold),OPTIONS,P,M,U,h_bold);
    end    
end


%%

flag_y=0;
% output
%==========================================================================
for i = 1:ns  % time points
 
    % output - g_vsdi(x)
    if isfield(M.hg,'vsdi')      
        y(:,i) = spm_vec(g(spm_unvec(x(i,:),M.x),U.u(i,:),P,M));
        flag_y=1;
        
    end
    
    % output - g_cai(h.cai)
    if isfield(M.hg,'cai') %isfield(M.fH,'h_cai')      
        y_cai(:,i) = spm_vec(g_cai(spm_unvec(x_cai(i,:),M.hx.x_cai),U.u(i,:),P,M));
        xh.x_cai   = x_cai;
        flag_y=1;        
    end
    
    % output - g_bold(h.bold)
    if isfield(M.hg,'bold') %isfield(M.fH,'h_bold')  
        y_bold(:,i) = spm_vec(g_bold(spm_unvec(x_bold(i,:),M.hx.x_bold),U.u(i,:),P,M));
        xh.x_bold   = x_bold;
        flag_y=1;        
    end

    % output - neural state dynamics x
    if flag_y==0
        y(:,i) = spm_vec(spm_unvec(x(i,:),M.x));
    end
    
end
 
% transpose
%--------------------------------------------------------------------------
if isfield(M.hg,'vsdi') || (flag_y==0)
    y      = real(y');

    % hidden signal
    if (size(y,2) ~= M.mC.nobs.vsdi), disp('chk!-vsdi'); end    
    id = unique(M.mC.obs.vsdi);
    %if(length(id) ~= M.mC.nobs.vsdi) % deleted 2021-05-15       
    y = y(:,id(id>0));         
    %end  % deleted 2021-05-15
    if isfield(M.obs,'mask'), y_vsdi = y; end % added jykang 2021-01-26 for step control
        
else
    y      = [];    
end

if isfield(M.fH,'h_cai')
    y_cai   = real(y_cai');    

    % hidden signal
    if (size(y_cai,2) ~= M.mC.nobs.cai), disp('chk!-cai'); end    
    id = unique(M.mC.obs.cai);
    %if(length(id) ~= M.mC.nobs.cai)        % deleted 2021-05-15
    y_cai = y_cai(:,id(id>0));
    %end % deleted 2021-05-15
    
    y = [y, y_cai];
end

if isfield(M.fH,'h_bold')  %if ~isfield(M.fH,'h_cai'), M.fH.h_cai=[]; end 
    y_bold   = real(y_bold');
    
    % hidden signal
    if (size(y_bold,2) ~= M.mC.nobs.bold), disp('chk!-bold'); end    
    id = unique(M.mC.obs.bold);
    %if(length(id) ~= M.mC.nobs.bold) % deleted 2021-05-15
    y_bold = y_bold(:,id(id>0));
    %end    % deleted 2021-05-15
    
    y = [y, y_bold];
end

% added jykagn 2021-01-26
if isfield(M.obs,'mask')
    y_bf_mask = y;         
    y = [];
    if isfield(M.hg,'vsdi') 
        try % todo : delete this
            y_mask = sig_mask(M.obs.vsdi.dt, M.obs.dt, y_vsdi);
        catch
            y_mask = sig_mask(M.obs.vsdi.dt, M.dt, y_vsdi);
        end
        y = [y(:); y_mask(:)];
    end
    if isfield(M.hg,'cai')
        try % todo : delete this 
            y_mask = sig_mask(M.obs.cai.dt, M.obs.dt, y_cai);
        catch
            y_mask = sig_mask(M.obs.cai.dt, M.dt, y_cai);
        end
        y = [y(:); y_mask(:)];
    end  
    if isfield(M.hg,'bold')
        try % todo : delete this
            y_mask = sig_mask(M.obs.bold.dt, M.obs.dt, y_bold);
        catch
            y_mask = sig_mask(M.obs.bold.dt, M.dt, y_bold);
        end
        y = [y(:); y_mask(:)];
    end  
    
    %y = spm_vec(y);
end

xstat=x;

end
