clear variables;
close all;

tmin=0;tmax=1;
f=@(t,y)((t.^2).*exp(t)-y.*(3*t.^2-1));%fonction y'(t)
y0=-1;
yExact=@(t)(exp(t).*(1-4*exp(-t.^3))/3);%solution exact du problème
    

% solution approchï¿½e de l'eq. diff.
% 1. mï¿½thode d'Euler (h=0.1)
h=0.1;
[yEuler1,t1]=fct_Euler(y0,tmin,tmax,h,f);
eps1=abs(yEuler1-yExact(t1));   % erreur



% 2. mï¿½thode d'Euler (h=0.05)
h=0.05;
[yEuler2,t2]=fct_Euler(y0,tmin,tmax,h,f);
eps2=abs(yEuler2-yExact(t2));   % erreur


% 1. mï¿½thode RK2 (h=0.1 et beta=1)
h=0.05;beta=1;
[yRK,t3]=fct_RK2(y0,tmin,tmax,h,beta,f);
eps3=abs(yRK-yExact(t3));       % erreur



%gestion de l'affichage des courbes
figure(1);
subplot(1,2,1);
hold on;
plot(t1,yEuler1,'*-');
plot(t1,yExact(t1));
plot(t3,yRK,'*-');
plot(t2,yEuler2,'*-');
legend('Euler pas de 0.1','Solution Exacte','Runge-Kutta pas de 0.05','Euler pas de 0.05');
title('Solution exacte, Solution approchée par méthode d Euler et Rung-Kutta 2');
xlabel(' t ');
ylabel ('y(t)');

subplot(1,2,2)
hold on;
plot(t1,eps1,'x-');
plot(t2,eps2,'x-');
plot(t3,eps3,'x-');
legend('Courbe erreur Euler pas de 0.1','Courbe erreur Euler pas de 0.05','Courbe erreur Runge-Kutta pas de 0.05');
title('Courbes des erreurs associées aux solutions approchées');
xlabel(' t ');
ylabel ('erreur');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
maxerr = [];
hi = [];
for h = 0.01 : 0.01 : 0.1
    hi = [hi, h]; % tableau des différentes valeurs du pas
    [yRK4,t4]=fct_RK4(y0,tmin,tmax,h,f); %calcul de la solution
    maxerr = [maxerr, max(yRK4-yExact(t4))]; %calcul de l'erreur maximale
end

% affichage de l'erreur en fonction du pas
figure(2);hold on;
plot(hi, maxerr);
title('courbe de l erreur en fonction du pas avec méthode Runge-Kutta ordre 4');
legend('erreur maximale entre la solution et la solution approchée');
ylabel('erreur');
xlabel(' pas '); 

%%détermination du polynome à l'aide de polyfit 
somme=0;
imax=9;
rms=[];

for i=1:1:imax 
   P = polyfit(hi, maxerr, i);
   somme=(maxerr(i)-P(i)).^2;
   rms=[rms, sqrt(somme./imax)];
end

figure(3);
plot([1:imax],rms);
title('courbe de l ecart quadratique moyen en fonction de l ordre du polynome pour la méthode Runge-Kutta ordre 4');
legend('ecart quadratique moyen');
ylabel('valeur de l ecart quadratique moyen');
xlabel('orde du polynome '); 