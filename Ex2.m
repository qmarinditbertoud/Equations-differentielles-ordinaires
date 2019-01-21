% CPE - 4ETI - 28/02/2017
% Module "Analyse Num?rique"
% EDO-EDP
% Travaux pratiques
% Exercice 2 : r?solution de l'?q. diff. de Lokta-Volterra (mod?le
% proies-pr?dateurs : (x',y')=(f(t,x,y),g(t,x,y))) par les m?thodes d'Euler
% et de Runge-Kutta d'ordre 2, visualisation des r?sultats sous forme de 
% trajectoire tangente au champ de vecteurs defini par (x,y)->(f(t,x,y),g(t,x,y)))


clear;close all;

tmin=0;
tmax=15;

% param?tres du mod?le proies-pr?dateurs
alpha=1;  % taux de reproduction des proies
beta=0.5; % taux de mortalit? des proies
gamma=2;  % taux de reproduction des pr?dateurs
delta=1;  % taux de mortalit? des pr?dateurs

% conditions initiales
x0=2;     % proies
y0=0.5;   % predateurs

% second membre de l'?qu. diff. (x',y')=(f(t,x,y),g(t,x,y))
f=@(t,x,y)(x.*(alpha - beta*y));
g=@(t,x,y)(-y.*(gamma - delta*x));

% Calcul des populations des proies et des pr?dateurs
h=0.01;   % pas temporel

% 1. méthode d'Euler
[xEuler,yEuler,t]=Euler_2D(x0,y0,tmin,tmax,f,g,h);

% 2. méthode RK2
beta=1;
[xRK2,yRK2,t1]=RK2_2D(x0,y0,tmin,tmax,h,beta,f,g);
%modification des Conditions initiales
[xRK2bis,yRK2bis,t1]=RK2_2D(3,2,tmin,tmax,h,beta,f,g);



% Affichage des populations des proies et des pr?dateurs
% en fonction du temps
figure(1);hold on;
plot(t,xEuler,'r');
plot(t,yEuler,'b');
plot(t1,xRK2,'g');
plot(t1, yRK2, 'm');


title('affichage des populations des proies et des prédateurs');
xlabel('base de temps');
ylabel('valeur des populations');
legend('affichage proies par la méthode Euler ','affichage prédateurs par la méthode Euler','affichage proies par la methode Runge-Kutta ordre2 x0=2 ', 'affichage prédateurs par la méthode Runge-Kutta ordre2 y0=0.5 ');

figure(2);hold on;
plot(t1,xRK2,'g');
plot(t1, yRK2, 'm');
plot(t1,xRK2bis,'y');
plot(t1, yRK2bis, 'k');
title('affichage des populations des proies et des prédateurs');
xlabel('base de temps');
ylabel('valeur des populations');
legend('affichage proies par la methode Runge-Kutta ordre2 x0=2 ', 'affichage prédateurs par la méthode Runge-Kutta ordre2 y0=0.5 ','affichage proies par la methode Runge-Kutta ordre2 x0=2 ', 'affichage prédateurs par la méthode Runge-Kutta ordre2 y0=3');


% affichage de la trajectoire proies-pr?dateurs (tangente au champ de
% vecteurs d?fini par la fonction (x,y)->(f(t,x,y),g(t,x,y))
figure(3);hold on;

% champ de vecteurs (x,y)->(f(t,x,y),g(t,x,y))
N=40;
ux=linspace(0,7,N);
uy=linspace(0,7,N);
[x,y]=meshgrid(ux,uy);       % grille de coordonn?es (ux,uy)
fxy=f(t,x,y);gxy=g(t,x,y);   % calcul du champ de vecteurs
norme=(fxy.^2+gxy.^2).^0.5;  % normalisation des vecteurs
fxy=fxy./norme;gxy=gxy./norme;
quiver(x,y,fxy,gxy);         % affichage de fxy et gxy sous forme
                             % de champ de vecteurs


% trajectoires proies-pr?dateurs
% 1. m?thode de d'Euler

plot(xEuler,yEuler, 'r');

% 2a. méthode RK2 

%plot(xRK2, yRK2, 'g');

plot(x0,y0,'r*');

title('trajectoires proies-prédateurs avec x0=2 et y0=0.5');
xlabel('nombre de proies');
ylabel('nombre de prédateurs');
legend('champ de vecteurs','trajectoires pour la méthode d Euler');


% 2b. méthode RK2 et nouvelles conditions initiales
[xRK2bis,yRK2bis,t1]=RK2_2D(3,2,tmin,tmax,h,beta,f,g);

figure(4);hold on;
N=40;
ux=linspace(0,7,N);
uy=linspace(0,7,N);
[x,y]=meshgrid(ux,uy);       % grille de coordonn?es (ux,uy)
fxy=f(t,x,y);gxy=g(t,x,y);   % calcul du champ de vecteurs
norme=(fxy.^2+gxy.^2).^0.5;  % normalisation des vecteurs
fxy=fxy./norme;gxy=gxy./norme;
quiver(x,y,fxy,gxy);         % affichage de fxy et gxy sous forme
                             % de champ de vecteurs
plot(xRK2, yRK2, 'g');
plot(xRK2bis, yRK2bis, 'r');
                             
plot(x0,y0,'g*');
plot(3,2,'r*');

title('trajectoires proies-prédateurs ');
xlabel('nombre de proies');
ylabel('nombre de prédateurs');
legend('champ de vecteurs','trajectoires pour la méthode Runge-Kutta pour x0=2 et y0=0.5','trajectoires pour la méthode Runge-Kutta pour x0=3 et y0=2');

