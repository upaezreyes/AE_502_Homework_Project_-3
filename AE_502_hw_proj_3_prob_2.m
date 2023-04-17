% AE 402: Orbital Mechanics
% Homework #2
% Newton's Method
function AE_502_hw_proj_3_prob_2

clc
clear variables 

pi = 3.1416; % pi
mu = 1; % AU^3/TU^2, gravitational constant
a = 1; % AU
n = sqrt(mu./(a.^3)); 
i = 45.*(pi./180); % inclination angle 
e = 0.5; % eccentricity 

% Given Orbital Elemens: 
L0 = sqrt(mu.*a); % L initial
G0 = L0.*sqrt(1 - (e.^2)); % G initial
H0 = G0.*cos(i); % H initial 
g0 = 0; % g = w, initial 
MA0 = 0; % l = M, mean anomaly initial 
h0 = 0; % h = omega, initial 

T0 = 2.*pi.*sqrt((a.^3)./mu); % sec, initial period 

% store initial orbital elements: 
coe0 = [L0 G0 H0 g0 MA0 h0]; 


% Use ODE45 to integrate the Gauss variational equations (Equations
% 12.89) from t0 to tf:
t0 = 0; % sec, initial time
days = 100; % days of time evolution
tf = days;  % sec, final time
nout = 1000; %number of solution points to output for plotting purposes
tspan = linspace(t0, tf, nout); % time steps

options = odeset(...
 'reltol', 1.e-8, ... 
 'abstol', 1.e-8, ...
'initialstep', T0/1000);
y0 = coe0';
[t,y] = ode45(@rates, tspan, y0, options);

% Assign the time histories mnemonic variable names:
L = y(:,1);
G = y(:,2);
H = y(:,3);
g = y(:,4);
MA = y(:,5);
h = y(:,6);

%a = (h.^2)./(mu.*(1-e.^2)); 
%disp(l)

%Plot the time histories of the osculating elements:
figure
subplot(3,1,1)
plot(t,h*(180./pi))
title('Right Ascension (degrees)')
xlabel('hours')
ylabel('degrees')
grid on
subplot(3,1,2)
plot(t,g*(180./pi))
title('Argument of Perigee (degrees)')
xlabel('hours')
ylabel('degrees')
grid on
subplot(3,1,3)
plot(t,MA)
title('Mean Anomaly')
xlabel('hours')
ylabel('')
grid on

figure
subplot(3,1,1)
plot(t,L)
title('L')
xlabel('hours')
ylabel('')
grid on
subplot(3,1,2)
plot(t,G)
title('G')
xlabel('hours')
ylabel('')
grid on
subplot(3,1,3)
plot(t,H)
title('H')
xlabel('hours')
ylabel('')
grid on
% ==========================================
% plotting r(t) = x + y + z 
% ===========================================
t = linspace(0,100,1000); % time units 
M  = sqrt(mu./(a.^3)).*t; % Mean anomaly angle

%M = linspace(0,-50.6124,1000); 
%M = MA'; 
Eo = M; % approximation initial value of Eo = M 


nmax = 100; % number of iterations
E = Eo; % initial value of E 

F =  @(Eo) M - Eo + e.*sin(Eo); % Function: M - E + e*sin(E)
Fdx = @(Eo) - 1 + e.*cos(Eo);  % dF/dE: -1 + e*cos(E)

% K = t./T0; 
% F =  @(Eo) M - (E - Eo) + e.*(sin(E) - sin(Eo)) - 2.*pi.*K; % Function: M - E + e*sin(E)
% Fdx = @(Eo) -1 + e.*cos(E);  % dF/dE: -1 + e*cos(E)

Evals = Eo; % array of iterations

for i=1:nmax
    E = E - F(E)/Fdx(E); % Newton's Method formula
    %Evals = [Evals;E]; % Values of E in each iterations
    
end

omega = h ;  
w_p = g; 
A = sqrt((1+e)./(1-e)).*tan(E./2); 
f = 2.*atan(A); % radians, true anomaly 
theta = w_p + f; % radians, argument of latitude at epoch
r = (a.*(1-(e.^2)))./(1 + e.*cos(f)); % km, distance 

r_x = r.*(cos(omega).*cos(theta) - sin(omega).*sin(theta).*cos(i)); 
r_y = r.*(sin(omega).*cos(theta) + cos(omega).*sin(theta).*cos(i)); 
r_z = r.*(sin(theta).*sin(i)); 

r_t2 = [r_x; r_y; r_z]; 
r_t2_mag = sqrt(sum(r_t2.^2)); 

% figure
% plot(t, r_t2)

% figure
% plot(t, r_t2_mag)

% figure
% plot3(r_x, r_y, t)
% xlabel("x-axis (AU)")
% ylabel("y-axis (AU)")
% zlabel("Time (TU)")
% title("Position r(t) vs Time (TU)")

figure
plot3(r_x, r_y, r_z)
xlabel("x-axis (AU)")
ylabel("y-axis (AU)")
zlabel("z-axis (AU)")
title("Position r(t) vs Time (TU)")

%=================================================================
% Plotting the orbital by using the Delaunay varibles 
%==============================================================

x = r.*cos(f + g).*cos(h) - (H./G).*r.*sin(f + g).*sin(h); 
y = r.*cos(f + g).*sin(h) + (H./G).*r.*sin(f + g).*cos(h); 
z = (sqrt(G.^2 - H.^2)./G).*r.*sin(f + g); 

figure
plot3(x, y, z)
xlabel("x-axis (AU)")
ylabel("y-axis (AU)")
zlabel("z-axis (AU)")
title("Position r(t) vs Time (TU)")


%============================================================
% f and y functions: 
%==============================================================

% initial position, r0: 
r0_v = [a.*(cos(Eo) - e); a.*sqrt(1-e.^2).*sin(Eo)];
r0 = sqrt(sum(r0_v.^2)); 

% initial velocity, v0: 
v0_v = [-(sqrt(mu.*a)./r0).*sin(Eo); (sqrt(mu.*a)./r0).*sqrt(1-e.^2).*cos(Eo)]; 
v0 = sqrt(sum(v0_v.^2)); 

% f & g functions: 
f_f = (a./r0).*((cos(E)-e).*cos(Eo) + sin(E).*sin(Eo)); 
g_f = sqrt((a.^3)./mu).*(sin(E-Eo) - e.*(sin(E)-sin(Eo)));

% f_f = 1 - (a./r).*(1 - cos(E - Eo));
% g_f = t - sqrt((a.^3)./mu).*(sin(E-Eo) - e.*(sin(E)-sin(Eo)));

% position, r(t): 
r_t = f_f.*r0_v + g_f.*v0_v; 
r_t_mag = sqrt(sum(r_t.^2)); 
%disp(r_t)


% plots: 

r_ta = f_f.*a.*(cos(Eo) - e) - g_f.*(sqrt(mu.*a)./r0).*sin(Eo); 
r_tb = f_f.*a.*sqrt(1-e.^2).*sin(Eo) + g_f.*(sqrt(mu.*a)./r0).*sqrt(1-e.^2).*cos(Eo);

figure 
plot3(r_ta,r_tb, t)
xlabel("x-axis (AU)")
ylabel("y-axis (AU)")
zlabel("time (t)")
title("Position, r(t) vs Time (t), f & g functions")

% Subfunction:
%====================================================================
function dfdt = rates(~,f)
%w===================================================================
%
% This function calculates the time rates of the orbital elements
% from Gauss’s variational equations (Equations 12.89).
%-------------------------------------------------------------------------
% The orbital elements at time t:
L = f(1);
G = f(2);
H = f(3);
g = f(4);
MA = f(5);
h = f(6);

% r = h^2/mu/(1 + e*cos(M)); % the radius
% u = w + M; % argument of latitude
w = 0.01; % ratating speed

% Delaunay element rates: 
Ldot = 0; 
Gdot = 0; 
Hdot = 0; 
gdot = 0; 
%MAdot = -1./(2.*L.^2); 
MAdot =  1./(L.^3); 
hdot = w; 

% pass these rates back to ODE45 in the array dfdt:
dfdt = [Ldot Gdot Hdot gdot MAdot hdot]';
end

end 
