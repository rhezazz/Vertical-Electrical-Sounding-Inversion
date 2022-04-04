clear 
clc
format long
res = [50 500 100];
thk = [5 15];
n = 0:17;
AB2 = 10.^(n/6);
rho_app = VES1D(res,thk,AB2);
rho_app = rho_app + rho_app.*normrnd(0,0.1,[1,length(rho_app)]);
d = [AB2' rho_app'];
writematrix(d,"Kurva K Noise 10%.dat",'delimiter','tab')
