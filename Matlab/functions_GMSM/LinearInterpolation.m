% Linear Interpolation of Excitation Numerical SDF Solution Method
% Function to calculate displacement & velocity response histories
%
% Written by: Neal (Simon) Kwong
% Date Modified: 04 Oct 2012
%% ----------------------DEFINITIONS---------------------------------------
% flag = indicator for ground motion (1 = yes)
% wn   = natural circular frequency in rad/sec^2
% zeta = damping ratio
% u0   = initial displacement
% v0   = initial velocity
% p    = excitation vector/file (either p or uddg - not peff)
% dt   = time step in excitation file
% m    = mass of SDF system
% wn2  = square of wn
% k    = stiffness of SDF system
% N    = number of data points in excitation file
% u    = relative displacement history
% v    = relative velocity history
%
% Note: Units are specified by excitation file p
%% ------------------------------------------------------------------------
function [u,v] = LinearInterpolation(flag,wn,zeta,u0,v0,p,dt,k)
%% Computation of constants for procedure
wn2 = wn^2;
m = k/wn2;
if flag == 1
    p = -m*p; % Convert p=uddg to be peff
end
wd = wn*sqrt(1-zeta^2);
N = length(p);
A = exp(-zeta*wn*dt)*((zeta/sqrt(1-zeta^2))*sin(wd*dt)+cos(wd*dt));
B = exp(-zeta*wn*dt)*sin(wd*dt)/wd;
C = 1/k*(2*zeta/(wn*dt)+exp(-zeta*wn*dt)*(((1-2*zeta^2)/(wd*dt)-zeta/sqrt(1-zeta^2))*sin(wd*dt)-(1+2*zeta/(wn*dt))*cos(wd*dt)));
D = 1/k*(1-2*zeta/(wn*dt)+exp(-zeta*wn*dt)*((2*zeta^2-1)/(wd*dt)*sin(wd*dt)+2*zeta/(wn*dt)*cos(wd*dt)));
Ap = -exp(-zeta*wn*dt)*(wn/sqrt(1-zeta^2)*sin(wd*dt));
Bp = exp(-zeta*wn*dt)*(cos(wd*dt)-(zeta/sqrt(1-zeta^2))*sin(wd*dt));
Cp = 1/k*(-1/dt+exp(-zeta*wn*dt)*((wn/sqrt(1-zeta^2)+zeta/(dt*sqrt(1-zeta^2)))*sin(wd*dt)+cos(wd*dt)/dt));
Dp = 1/(k*dt)*(1-exp(-zeta*wn*dt)*(zeta/sqrt(1-zeta^2)*sin(wd*dt)+cos(wd*dt)));

%% Initialization of response histories
u = zeros(1,N);     v = zeros(1,N);
u(1) = u0;          v(1) = v0;

% Superposition of Exact Solutions as Recurrence Formulas
for i=1:(N-1)
    u(i+1) = A*u(i)+B*v(i)+C*p(i)+D*p(i+1);
    v(i+1) = Ap*u(i)+Bp*v(i)+Cp*p(i)+Dp*p(i+1);
end

end
