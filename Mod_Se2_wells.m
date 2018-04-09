% Se2_wells.m.; % 3 wells to form a quasi -FET
clear all
NN = 250; % Number of points in the problem space.
hbar = 1.054e-34; % Plank ?s constant
melec = 9.1e-31; % Mass of an electron
eV2J = 1.6e-19; % Energy conversion factors
J2eV = 1/eV2J;
del_x = .2e-9; % The cells size
dt = 1e-16; % Time steps
ra = (0.5 *hbar/melec)*(dt/del_x^2) % ra must be < .1
DX = del_x *1e9; % Cell size in nm.
XX = (DX:DX:DX *NN); % Length in nm for plotting
% --- Specify the coupled wells -----------------------------
V = zeros(1,NN);
for n =99:101
V(n) = .7 *eV2J;
end
for n =126:128
V(n) = .7 *eV2J;
end
% Middle well
for n =102:125
V(n) = -0.005*eV2J; % Mimics a gate voltage
end
% ----------- Initialized a waveform ------------------------
lambda = 200; % Pulse wavelength
prl = zeros(1,NN);
pim = zeros(1,NN);
ptot = 0.;
for n =2:99
prl(n) = sin(4 *pi*n/lambda) ;
ptot = ptot + prl(n)^2 + pim(n) ^2;
end
pnorm = sqrt(ptot);
% Normalize and check
ptot = 0.;
for n =1:NN
prl(n) = prl(n)/pnorm;
pim(n) = pim(n)/pnorm;
ptot = ptot + prl(n)^2 + pim(n) ^2;
end
ptot
T = 0;
n_step = 1;
while n_step > 0
n_step = input( 'How many time steps -->');
% ------------ This is the main FDTD loop ------------
for m =1:n_step
T = T + 1;
for n =2:NN-1
prl(n) = prl(n) - ra *(pim(n-1) -2*pim(n) + pim(n +1)) ...
+ (dt/hbar) *V(n)*pim(n);
end
for n =2:NN-1
pim(n) = pim(n) + ra *(prl(n-1) -2*prl(n) + prl(n +1)) ...
- (dt/hbar) *V(n)*prl(n);
end
end
% ------------------------
% Check normalization
ptot = 0.;
for n =1:NN
ptot = ptot + prl(n)^2 + pim(n) ^2;
end
ptot
% Calculate the expected values
PE = 0.;
for n =1:NN
psi(n) = prl(n) + 1j *pim(n);
PE = PE + psi(n) *psi(n)'*V(n);
end
PE = PE *J2eV;
ke = 0. + 1j * 0.;
for n =2:NN-1
lap_p = psi(n +1) - 2 *psi(n) + psi(n -1);
ke = ke + lap_p *psi(n)';
end
KE = -J2eV*((hbar/del_x)^2/(2*melec))*real(ke);
%subplot(3,1,1)
plot(XX,prl,'k')
hold on
plot(XX,pim,'-.k')
plot(XX,J2eV*V,'--k')
hold off
axis( [ 0 DX *NN -0.2 .75 ])
TT = text(5,.15,sprintf( '%7.2f ps ',T*dt*1e12));
set(TT,'fontsize',12)
TT = text(35,.15,sprintf( 'KE = %5.2f meV ',1e3*KE));
set(TT,'fontsize',12)
xlabel('nm')
TT = ylabel( 'y','FontName','Symbol','fontsize',12);
set(gca,'fontsize',12)
title('Se2-wells')
T
saveas(gcf,'se2_wells.png')      % This saves the picture to a file
end