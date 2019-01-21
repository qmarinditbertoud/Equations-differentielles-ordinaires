function[y,t]=fct_RK2(y0,tmin,tmax,h,beta,f)

t=tmin:h:tmax;
y=zeros(1,length(t));
y(1)=y0;

for k=2 : length(t)
    H=h/(2*beta)
    y(k)=y(k-1)+h*( (1-beta) * f(t(k-1),y(k-1)) + beta* f((t(k-1)+H),y(k-1)+(H*f(t(k-1),y(k-1)))) );
    
end
end

  