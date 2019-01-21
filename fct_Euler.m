function[y,t]=fct_Euler(y0,tmin,tmax,h,f)

t=tmin:h:tmax;%création d'un tableau pour la variable t avec pas de h
y=zeros(1,length(t));%tableau de 0
y(1)=y0;%initialisation du premier terme d'indice 1 = à y0
for k=2 : length(t) %boucle calculant les termes d'indices superieur
    y(k)=y(k-1)+h*f(t(k-1),y(k-1));
end
end