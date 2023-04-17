function Bonus_question

clc
clear variables 

% a = linspace(1,10,10); 
% a_length = length(a); 

for j = 1:1:20 % this loop increases by 1 every iteration
               % j represents the semi-major axis

pi = 3.1416; % pi
mu = 1; % AU^3/TU^2, gravitational constant
%a = 1; % AU

n = sqrt(mu./(1.^3)); % mean motion
i =45.*(pi./180); % inclination angle 
e = 0.5; % eccentricity
w = 0.5; % rotating speed

% Given Orbital Elemens: 
L0 = n.*j.^2; % L initial
G0 = L0.*sqrt(1 - (e.^2)); % G initial
H0 = G0.*cos(i); % H initial 
g0 = 0; % g = w, initial 
MA0 = 0; % l = M, mean anomaly initial 
h0 = 0; % h = omega, initial 

T0 = 2.*pi.*sqrt((j.^3)./mu); % sec, initial period

% store initial orbital elements: 
coe0 = [L0 G0 H0 g0 MA0 h0 L0 G0 H0 g0 MA0 h0]; 


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
% Delaunay variables of the old Hamiltonina used in question 2: 
L = y(:,1);
G = y(:,2);
H = y(:,3);
g = y(:,4);
MA = y(:,5);
h = y(:,6);

% put the values into a vector 
L_v(:,j) = L; 
G_v(:,j) = G; 
H_v(:,j) = H; 
g_v(:,j) = g; 
MA_v(:,j) = MA; 
h_v(:,j) = h; 

% Delaunay variables of the new Hamiltonian afte the 
% Hori-Lie-Deprit was applied: 
L1 = y(:,7);
G1 = y(:,8);
H1 = y(:,9);
g1 = y(:,10);
MA1 = y(:,11);
h1 = y(:,12);

% put the values into a vector 
L1_v(:,j) = L1; 
G1_v(:,j) = G1; 
H1_v(:,j) = H1; 
g1_v(:,j) = g1; 
MA1_v(:,j) = MA1; 
h1_v(:,j) = h1; 

% new Delaunay variables (L',G',H',...): 
L1_p = L1 + w.*(L1.^3).*H1; 
G1_p = G1; 
H1_p = H1; 
MA1_p = MA1 - 3.*w.*(L1.^2).*H1.*MA1; 
g1_p = g1; 
h1_p = h1 - w.*(L1.^3).*MA1; 


%==========================================================
% equinoctial elements: h, k, p, and q 
%=====================================================

% h and k: 
he1 = e.*sin(g1_p + h1_p); 
he2 = e.*sin(g + h);
ke1 = e.*cos(g1_p + h1_p); 
ke2 = e.*cos(g + h); 

% put the values into vectors 
he1_v(:,j) = he1; 
he2_v(:,j) = he2; 
ke1_v(:,j) = ke1; 
ke2_v(:,j) = ke2;

% p and q 
pe1 = tan(i./2).*sin(h1_p); 
pe2 = tan(i./2).*sin(h); 
qe1 = tan(i./2).*cos(h1_p); 
qe2 = tan(i./2).*cos(h); 

% put the values into vectors 
pe1_v(:,j) = pe1; 
pe2_v(:,j) = pe2; 
qe1_v(:,j) = qe1; 
qe2_v(:,j) = qe2;

end % end of the loop 

%================================================================
% Plotting Delaunay variablies: 
%=============================================================
% figure
% plot(t,h_v(:,1)*(180./pi))
% for j = 2:1:20
%     hold on
%     plot(t,h_v(:,j))
% end
% hold on 
% plot(t,h1_v(:,1)*(180./pi),'--')
% for j = 2:1:20
%     hold on
%     plot(t,h1_v(:,j),'--')
% end
% title('Right Ascension (degrees)')
% xlabel('hours')
% ylabel('degrees')
% legend("Q2","Q1")
% grid on
% 
% figure
% plot(t,g_v(:,1)*(180./pi))
% for j = 2:1:20
%     hold on
%     plot(t,g_v(:,j))
% end
% hold on 
% plot(t,g1_v(:,1)*(190./pi),'--')
% for j = 2:1:20
%     hold on
%     plot(t,g1_v(:,j),'--')
% end
% title('Argument of Perigee (degrees)')
% xlabel('hours')
% ylabel('degrees')
% legend("Q2","Q1")
% grid on
% 
% figure
% plot(t,MA_v(:,1)*(180./pi))
% for j = 2:1:20
%     hold on
%     plot(t,MA_v(:,j))
% end
% hold on 
% plot(t,MA1_v(:,1)*(180./pi),'--')
% for j = 2:1:20
%     hold on
%     plot(t,MA1_v(:,j),'--')
% end
% title('Mean Anomaly')
% xlabel('hours')
% ylabel('')
% legend("Q2","Q1")
% grid on
% 
% figure
% plot(t,L_v(:,1))
% for j = 2:1:20
%     hold on
%     plot(t,L_v(:,j))
% end
% hold on 
% plot(t,L1_v,'--')
% for j = 2:1:20
%     hold on
%     plot(t,L1_v(:,j),'--')
% end
% title('L')
% xlabel('hours')
% ylabel('')
% legend("Q2","Q1")
% grid on
% 
% figure
% plot(t,G_v(:,1)*(180./pi))
% for j = 2:1:20
%     hold on
%     plot(t,G_v(:,j))
% end
% hold on 
% plot(t,G1_v(:,1)*(180./pi),'--')
% for j = 2:1:20
%     hold on
%     plot(t,G1_v(:,j),'--')
% end
% title('G')
% xlabel('hours')
% ylabel('')
% legend("Q2","Q1")
% grid on
% 
% figure
% plot(t,H_v(:,1)*(180./pi))
% for j = 2:1:20
%     hold on
%     plot(t,H_v(:,j))
% end
% hold on 
% plot(t,H1_v(:,1)*(180./pi),'--')
% for j = 2:1:20
%     hold on
%     plot(t,H1_v(:,j),'--')
% end
% title('H')
% xlabel('hours')
% ylabel('')
% legend("Q2","Q1")
% grid on

%=======================================================
% Plotting equinoctial elements : h vs q and p vs q
%=======================================================

% h vs k: 
figure 
plot(he1_v(:,1),ke1_v(:,1),LineWidth = 3)
for j = 2:1:20
    hold on
    plot(he1_v(:,j),ke1_v(:,j),'--',LineWidth = 3)
end
hold on 
plot(he2_v(:,1),ke2_v(:,1),LineWidth = 3)
for j = 2:1:20
    hold on
    plot(he2_v(:,j),ke2_v(:,j),'--',LineWidth = 3)
end
xlabel('h')
ylabel('k')
title('h vs k')

% p vs q: 
figure 
plot(pe1_v(:,1),qe1_v(:,1),LineWidth = 3)
for j = 2:1:20
    hold on
    plot(pe1_v(:,j),qe1_v(:,j),'--',LineWidth = 3)
end
hold on 
plot(pe2_v(:,1),qe2_v(:,1),LineWidth = 3)
for j = 2:1:20
    hold on
    plot(pe2_v(:,j),qe2_v(:,j),'--',LineWidth = 3)
end
xlabel('p')
ylabel('q')
title('p vs q')

% Subfunction:
%====================================================================
function dfdt = rates(~,f)
%w===================================================================
%
% This function calculates the time rates of the orbital elements
% from Gauss’s variational equations (Equations 12.89).
%-------------------------------------------------------------------------
% The orbital elements at time t:
% Delaunay variables of the old Hamiltonina: 
L = f(1);
G = f(2);
H = f(3);
g = f(4);
MA = f(5);
h = f(6);

Ldot = 0; 
Gdot = 0; 
Hdot = 0; 
gdot = 0; 
%MAdot = -1./(2.*L.^2); 
MAdot =  1./(L.^3); 
hdot = w; 

% Delaunay variables of the new Hamiltonian
L1 = f(7);
G1 = f(8);
H1 = f(9);
g1 = f(10);
MA1 = f(11);
h1 = f(12);

L1dot = 0; 
G1dot = 0; 
H1dot = 0; 
g1dot = 0; 
%MAdot = -1./(2.*L.^2); 
MA1dot = 1./(L.^3); 
h1dot = 0; 

% pass these rates back to ODE45 in the array dfdt:
dfdt = [Ldot Gdot Hdot gdot MAdot hdot L1dot G1dot H1dot g1dot MA1dot h1dot]';
end

end 

