
% CPE - 4ETI - 28/02/2017
% Module "Analyse Num�rique"
% EDO-EDP
% Travaux pratiques
% Exercice 2 : r�solution de l'�q. diff. de Lokta-Volterra (mod�le
% proies-pr�dateurs : (x',y')=(f(t,x,y),g(t,x,y))) par les m�thodes d'Euler
% et de Runge-Kutta d'ordre 2, visualisation des r�sultats sous forme de 
% trajectoire tangente au champ de vecteurs defini par (x,y)->(f(t,x,y),g(t,x,y)))


clear;close all;

tmin=0;
tmax=15;

% param�tres du mod�le proies-pr�dateurs
alpha=1;  % taux de reproduction des proies
beta=0.5; % taux de mortalit� des proies
gamma=2;  % taux de reproduction des pr�dateurs
delta=1;  % taux de mortalit� des pr�dateurs

% conditions initiales
x0=2;     % proies
y0=0.5;   % predateurs

% second membre de l'�qu. diff. (x',y')=(f(t,x,y),g(t,x,y))
f=@(t,x,y)((x.*(alpha-beta*y)));
g=@(t,x,y)((-y.*(gamma-delta*x)));

% Calcul des populations des proies et des pr�dateurs
h=0.01;   % pas temporel

% 1. m�thode d'Euler
[xEuler,yEuler,t]=Euler_2D(y0,x0,tmin,tmax,f,g,h);

% 2. m�thode RK2
%beta=1;
[xRK2,yRK2,t1]=RK2_2D(x0,y0,tmin,tmax,h,beta,f,g);


% Affichage des populations des proies et des pr�dateurs
% en fonction du temps
figure(1);hold on;

plot(t,xEuler,'*-');
plot(t,yEuler,'*-');




% affichage de la trajectoire proies-pr�dateurs (tangente au champ de
% vecteurs d�fini par la fonction (x,y)->(f(t,x,y),g(t,x,y))
figure(2);hold on;

%champ de vecteurs (x,y)->(f(t,x,y),g(t,x,y))
N=40;
ux=linspace(0,7,N);
uy=linspace(0,7,N);
[x,y]=meshgrid(ux,uy);       % grille de coordonn�es (ux,uy)
fxy=f(t,x,y);                % calcul du champ de vecteurs
gxy=g(t,x,y);   
norme=(fxy.^2+gxy.^2).^0.5;  % normalisation des vecteurs
fxy=fxy./norme;
gxy=gxy./norme;
quiver(x,y,fxy,gxy);         % affichage de fxy et gxy sous forme
                             % de champ de vecteurs
                             


% trajectoires proies-pr�dateurs
% 1. m�thode de d'Euler

plot(xEuler,yEuler);

% 2a. m�thode RK2 
figure(3);hold on;

plot(t1,xRK2,'*-');
plot(t1,yRK2,'*-');



% 2b. m�thode RK2 et nouvelles conditions initiales







% M�thode d'Euler
% function [x,y,t]=Euler_2D(x0,y0,tmin,tmax,pas,F,G)
%     x=0;
%     y=0;
%     t=0;
% 
% % M�thode de Runge-Kutta d'ordre 2
% function [x,y,t]=RK2_2D(x0,y0,tmin,tmax,pas,beta,F,G)
%     x=0;
%     y=0;
%     t=0;
%   