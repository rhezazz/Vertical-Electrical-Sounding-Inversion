%Very Fast Simulated Annealing for VES Data Inversion
%Mohammad Rheza Zamani
%Reference : Grandis,H.(2009): Pengantar pemodelan inversi geofisika, Jakarta:HAGI
%Reference : W. Srigutomo, M. Heriyanto, and M. Hilmi Aufa. Gravity Inversion of Talwani Model using Very Fast Simulated Annealing. Journal of Mathematical and Fundamental Sciences, Vol. 51, No. 2, 2019, 177-190. doi: 10.5614/j.math.fund.sci.2019.51.2.7.
tic;
clear all;
clc;
%Synthetic Model
k = 0:17;
AB2 = 10.^(k/6);
res = [500 250 1000];
thk = [5 25];
rho_app =VES1D(res,thk,AB2);
%Inversion Parameter
nlayer = 3; 
nitr = 500;
T = 1;
dec = 1;
%Generate random model
LBR = [1 1 1];
UBR = [2000 2000 2000];
LBT = [1 1];
UBT = [50 50];
rho1(1 , :) = [1000 1000 1000];
thick1(1, :) = [20 20 50];

%Calculate respon and misfit for each model
[rho_app_1] = VES1D(rho1(1,:),thick1(1,:),AB2);
app_rho1(1,:) = rho_app_1;
[misfit1] = misfit_VES(rho_app,app_rho1(1,:));
E1 = misfit1;

%Inversion process
for itr = 1 : nitr
    rho_int(1 , :) = LBR + rand*(UBR - LBR);
    thick_int(1, :) = LBT + rand*(UBT - LBT);
    ui = rand;
    yi = sign(ui-0.5)*T*((((1 + (1/T)))^abs(2*ui-1))-1);
    rho2(1 , :) = rho_int + yi*(UBR - LBR);
    thick2(1, :) = thick_int + yi*(UBT - LBT);
    [rho_app_2] = VES1D(rho2(1,:),thick2(1,:),AB2);
    app_rho2(1,:) = rho_app_2;
    [misfit2] = misfit_VES(rho_app,app_rho2(1,:));
    E2 = misfit2;

    delta_E = E2 -E1;
    if delta_E < 0
        rho1 = rho2;
        thick1 = thick2;
         E1 = E2;
    else
        P = exp((-delta_E)/T);
        if  P >= rand
           rho1 = rho2;
           thick1 = thick2;
           E1 = E2;
        end
    end
    [rho_app_new] = VES1D(rho1,thick1,AB2);
    Egen(itr)=E1;
    T = T*exp(-dec*(itr)^(1/(2*nlayer)-1));
    Temperature(itr) = T;
end

%Data vizualization 
r_plot = [0, res];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho1];
Depth_mod = [0,cumsum(thick1),max(thick1)*100];

figure(1)
subplot(1,6,[1 3])
loglog(AB2,rho_app_new,'r',AB2,rho_app,'ob','MarkerSize',6,'LineWidth',2.5);
axis([1 10^3 1 10^4]);
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('AB/2 (m)','FontSize',8,'FontWeight','Bold');
ylabel('App. Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}Respon  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
grid on
subplot(1,6,[5 6])
stairs(r_plot,t_plot,'--r','Linewidth',2.5);
hold on
stairs(r_mod,Depth_mod,'-b','Linewidth',1.5);
hold off
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','Southwest');
axis([1 10^4 0 50]);
xlabel('Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
set(gca,'YDir','Reverse');
set(gca, 'XScale', 'log');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RMSE','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

figure(3)
plot(1:nitr,Temperature,'b','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('Temperature','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Penurunan Temperature ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on
time = toc;