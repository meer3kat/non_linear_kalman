clear all

load('SP500');
%para = [u,k,th,lam,ro]
%       u, kappa, theta, lamda, ro
u = 0.2016;
k = 2.1924;
th = 0.0133;
lam = 0.294;
ro = -0.6143;
delta_t = 1/252;


x0 = 0.0233;
p0 = 0.00001;
z = log(data);

Fk = 1 - k*(delta_t) + 0.5*ro*lam*(delta_t);
Ak = k*th*delta_t - lam*ro*u*delta_t;
Qk = lam*lam*(1 - ro*ro)*x0*delta_t;

Hk = -0.5*delta_t;
Bk = z(1) + u*delta_t;
Rk = x0 * delta_t;

xp(1) = Fk*x0 + Ak; %first prediction of x
%zp(1) = Hk*xp(1) + Bk; %predicted z (log price)
pp(1) = Fk*p0*Fk'+Qk;


Kk = pp(1)*Hk'*(Hk*pp(1)*Hk' + Rk)^-1;
zk(1) = Hk*xp(1)+Bk;
xm(1) = xp(1) + Kk*(z(1) - zk(1));
pm(1) = (1 - Kk*Hk)*pp(1);


for i = 2:1:length(data)
    
    Fk = 1 - k*(delta_t) + 0.5*ro*lam*(delta_t);
    Ak = k*th*delta_t - lam*ro*u*delta_t + lam*ro*(z(i)-z(i-1));
    Qk = lam * lam * (1-ro*ro)*xp(i-1)*delta_t;
    Hk = -0.5*delta_t;
    Bk = z(i-1) + u*delta_t;
    Rk = xp(i-1) * delta_t;
    
    xp(i) = Fk * xm(i-1) + Ak;
    pp(i) = Fk * pm(i-1) * Fk' + Qk;
    
    K(i) = pp(i)*Hk'*(Hk*pp(i)*Hk'+ Rk)^-1;
    
    zk(i) = Hk*xp(i) + Bk;%from volatility to price
    
    xm(i) = xp(i) + K(i)*(z(i) - zk(i));
    pm(i) = (1 - K(i)*Hk)*pp(i);
 
end

plot(dtime, z, 'b', dtime, zk, 'r-')
datetick('x','yyyy');