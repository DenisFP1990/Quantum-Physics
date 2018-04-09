% Programa para verificar que a una distancia del canal de 19.4 nm
% la partícula de 1keV atraviesa el MOSFET
% Se2_wells.m.
% 3 wells to form a quasi -FET
clear all
% Nuestro tamaño del MOSFET es ahora 59.4 pm aproximadamente
% ya que el tamaño del drain y source son 20 nm y el tamaño del canal
% es 19.4 nm según lo obtenido en el literal a.
% 59.4 pm = 297 posiciones del vector NN
NN = 297; % Number of points in the problem space.
hbar = 1.054e-34; % Plank ?s constant
melec = 9.1e-31; % Mass of an electron
eV2J = 1.6e-19; % Energy conversion factors
J2eV = 1/eV2J;
% Ajuste de espacio en escala de Pico Metros
% Cada vector corresponde 0.2 pm
del_x = 0.2e-12; % The cells size
% Ajuste de desplazamiento en el tiempo
% Cada paso de simulacion es de 10 yoctosegundos
dt = 10e-24; % Time steps
ra = (0.5 *hbar/melec)*(dt/del_x^2) % ra must be < 0.1
% Adecuamos la escala de la celda en picometros
DX = del_x *1e12; % Cell size in pm.
XX = (DX:DX:DX *NN); % Length in pm for plotting
% --- Specify the coupled wells -----------------------------
% La barrera de potencial inicial se mantiene en la posición 99, es decir
% a 19.8 pm del inicio, el ancho del canal ahora debera ser de 19.4 pm
% Según lo calculado, por tanto este ahora se desplazara 19.4/0.2 que
% Corresponde a 97 posiciones del vector, por tanto tendrá su presencia
% de 102 a 199 y despúes empezara la segunda barrera.
% Las barreras de potencial para ser consideradas como infinitas deben
% ser mucho mayores al valor de la energía de la partícula, por tanto si
% la energía de la partícula es 1keV, la barrera debería tener un valor
% significativamente mayor, utilizaremos 200 keV debido a que el programa
% original tendiendo 1meV utilizaba barreras de 200meV, 200 veces mas el
% valor de la energía de la partícula
V = zeros(1,NN);
for n =99:101
V(n) = 200000 *eV2J;
end
for n =200:202
V(n) = 200000 *eV2J;
end
% NO aplicamos potencial al canal
% Middle well
for n =102:199
V(n) = -0.000*eV2J; % Mimics a gate voltage
end
% ----------- Initialized a waveform ------------------------
% Mantenemos lamda de 400 correspondiente al primer estado cuántico
% de la partícula
lambda = 400; % Pulse wavelength
prl = zeros(1,NN);
pim = zeros(1,NN);
ptot = 0.;
% La partícula se inicializa en los primeros 19.8 pm
for n =2:99
prl(n) = sin(4 *pi*n/lambda);
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
% El programa tiene que adecuar las escalas de impresion nuevamente
% Debido a que todo se encuentra en pm, la impresión se tiene que dar en
% esa escala, y el avance del tiempo lo plantearemos en attosegundos
% debido a que se requiere de una cantidad considerable de tiempo en
% yoctosegundos para visualizar el efecto tunelamiento, para no tener
% valores muy altos en yocto o zepto segundos, utilizaremos attosegundos
axis( [ 0 DX *NN -.2 .25 ])
% los pasos temporales son de 10 yoctosegundos = 10e-24, para ajuste en
% atto segundos multiplicamos por 0.1e18, 1attosegundo = 1e-18 seg
TT = text(5,.15,sprintf( '%7.2f as ',T*dt*0.1e18));
set(TT,'fontsize',12)
TT = text(35,.15,sprintf( 'KE = %5.2f keV ',1e-3*KE));
set(TT,'fontsize',12)
xlabel('pm')
TT = ylabel( 'y','FontName','Symbol','fontsize',12);
set(gca,'fontsize',12)
title('Se2-wells')
T
saveas(gcf,'se2_wells.png')      % This saves the picture to a file
end