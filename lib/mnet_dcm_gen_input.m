function U=mnet_dcm_gen_input(type,t,dt,tstart,tend)
%% External input 
% type 1 : gaussian distribution 
% type 2 : square function
if type == 1
    [u]    = genU(t);
    U.u      = u / max(u);
    U.dt     = dt;
elseif type == 2
    u       = zeros(size(t))';
    u(tstart:tend) = 1;
    U.u     = u;
    U.dt    = dt;
elseif type == 3
    U.u = pinknoise(length(t),1);
    U.dt    = dt;
else
    disp('cound not find input type');
end
end


function [u]=genU(t)
if nargin<1,t=[1:600]; end
if nargin<2,tp=[1000];end

%tau=0.2;amp=10;
tau=5000;amp=10;
u=zeros(length(t),1); u2=u;
for i=1:length(t)
    u(i)=0;
    for j=1:length(tp)
        u(i)=u(i)+amp*exp(-((i-tp(j))^2)/tau);
    end
end
end 