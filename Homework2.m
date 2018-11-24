%%% ECE 550 Homework 2
%%% MOS capacitor electrostatics: n-substrate 
%%% Code created by Shaikh S. Ahmed // SIUC
%%% Edited by Steven Smith

%%% Clear workspace:
clc 
clear all
close all
%%% Inputs:
Tox = 2e-9;
Nd = 1e23;
q = 1.6e-19;
Ks = 11.7;
Kox = 3.95;
k = 1.38e-23;  %% J/K
T = 300;
eps0 = 8.85e-12;
Nc = 2e25;
Nv = 1e25;
Eg = 1.1;      %% in eV
Eg = 1.1*q     %% in Joules (J)
%%% Calculation:
ni = sqrt(Nc*Nv)*exp(-Eg/(2*k*T))  % Eg and kT are both in J
phiF = (k*T)*log(Nd/ni)    % in Joules                                  
Ldi = sqrt((k*T*Ks*eps0)/(2*q*q*ni)) %% Joule based calculation 
ii=1;
 % -ve gate voltage
for Vg = -2:0.01:0           %% Vg is potential (Volts) here
xaxis(ii) = Vg;
    for V = 1.2:-0.00005:-1.2 %% V is a numeric (iterative) variable for calculating Vs in Volts
lhs = Vg; %% LHS of the K.V.L.
u = -V*q; %% u is in Joules
Fs = sqrt(exp(phiF/k/T)*(u/k/T + exp(-u/k/T) - 1) - exp(-phiF/k/T)*(u/k/T -exp(u/k/T) + 1));
rhs = -(Ks*Tox)/Kox * (k*T)/q * 1/Ldi * Fs + V; %% RHS of the K.V.L. (-ve for -ve Gate voltage)
if(lhs/rhs >= 0.99 && lhs/rhs <= 1.01)
Vs = V; %% Vs (surface potential) is in Volts 
Fss = Fs;
end
end
Qs(ii) = (Ks*eps0)*(k*T)/q*1/Ldi*Fss/q; %% #/cm2;


yaxis(ii) = Vs;
ii = ii + 1;
end
 % +ve gate voltage
for Vg = 0:0.01:2
xaxis(ii) = Vg;
  for V = -1.2:0.00005:1.2
 %calculating Vs in Volts
  lhs = Vg;
  u = -V*q;
 %yaxis is surface potential (volts)
 % Vg is potential (Volts) here
 % V is a numeric (iterative) variable for
 % LHS of the K.V.L.
 % u is in Joules
Fs = sqrt(exp(phiF/k/T)*(u/k/T + exp(-u/k/T) - 1) - exp(-phiF/k/T)*(u/k/T - exp(u/k/T) + 1));
rhs = (Ks*Tox)/Kox * (k*T)/q * 1/Ldi * Fs + V; %% RHS of the K.V.L.(+ve for +ve Gate voltage)
if(lhs/rhs >= 0.99 && lhs/rhs <= 1.01)
Vs = V; %% Vs (surface potential) is in Volts 
Fss = Fs;
end
end
Qs(ii) = -(Ks*eps0)*(k*T)/q*1/Ldi*Fss/q; %% #/cm2
yaxis(ii) = Vs; % yaxis is surface potential (volts)
ii = ii + 1;
end
%%% Capacitance:

for n = 2:1:ii-2
    Vgate(n) = xaxis(n);
    cap(n) = -(Qs(n)-Qs(n-1))*q/(Vgate(n)-Vgate(n-1)); 
    if(cap(n)<1e-10)
        cap(n) = cap(n-1);
    end
end

%%%% Figure Plotting
figure('Color','White'); 
h1=semilogy(xaxis,abs(Qs),'r');
%of gate voltage
% Line options
set(h1,'linestyle','-');
set(h1,'linewidth',2)
% Figure options
set(gcf,'Colormap',pink);
% Axis options
axis tight;
set(gca,'color',[1 1 1]*0.9); 
set(gca,'fontsize',12);
set(gca,'layer','top');
set(gca,'linewidth',2);
xlabel('GATE VOLTAGE, V_G_S [V]','Fontsize',16); 
ylabel('|Q_S|, [#/m^3]','Fontsize',16);
% save
saveas(h1,'nMOScap-Qs-Vgs.fig');
print nMOScap-Qs-Vgs.eps -deps;

%%%% Figure Plotting
figure('Color','White');
h1=plot(xaxis,yaxis,'r');
% Line options
set(h1,'linestyle','-');
set(h1,'linewidth',2)
% Figure options
set(gcf,'Colormap',pink);
% Axis options
axis tight;
set(gca,'color',[1 1 1]*0.9);
set(gca,'fontsize',12);
% Plotting the absolute charge as a function
set(gca,'layer','top');
set(gca,'linewidth',2);
xlabel('GATE VOLTAGE, V_G_S [V]','Fontsize',16); 
ylabel('PHI_S, [V]','Fontsize',16);
% save
saveas(h1,'nMOScap-phiS-Vgs.fig');
print nMOScap-phiS-Vgs.eps -deps;

%%%% Figure Plotting
figure('Color','White');
h1=plot(Vgate,cap,'r');
% Line options
set(h1,'linestyle','-');
set(h1,'linewidth',2)
% Figure options
% Plotting the absolute charge
set(gcf,'Colormap',pink);
% Axis options
axis tight;
set(gca,'color',[1 1 1]*0.9); 
set(gca,'fontsize',12);
set(gca,'layer','top');
set(gca,'linewidth',2);
xlabel('GATE VOLTAGE, V_G_S [V]','Fontsize',16); 
ylabel('CAPACITANCE, [F/m^2]','Fontsize',16);
% save
saveas(h1,'nMOScap-Cap-Vgs.fig'); 
print nMOScap-Cap-Vgs.eps -deps;