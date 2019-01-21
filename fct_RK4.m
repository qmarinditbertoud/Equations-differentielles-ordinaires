% m?thode de Runge-Kutta d'ordre 4
function [y,t]=fct_RK4(y0,tmin,tmax,h,f)
  nbIters=floor((tmax-tmin)/h);
  y=zeros(1,nbIters+1);
  t=zeros(1,nbIters+1);
  y(1)=y0;
  t(1)=tmin;
  for k=1:nbIters
      k1=f(t(k),y(k));
      k2=f(t(k)+h/2,y(k)+(h/2)*k1);
      k3=f(t(k)+h/2,y(k)+(h/2)*k2);
      k4=f(t(k)+h,y(k)+h*k1);
      y(k+1)=y(k)+(h/6)*(k1+2*k2+2*k3+k4);
      t(k+1)=t(k)+h;
  end