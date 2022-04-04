%Simulated Annealing for VES Data Inversion
%Mohammad Rheza Zamani
%Reference : Kirkpatrick S, Gellat CD, Vecchi MP. Optimization by simulated annealing. Science 1983;220(4598):671-80.
clear all;
clc;
%Synthetic Model
n = 0:17;
AB2 = 10.^(n/6);
R = [20 100 10];
thk = [5 10];
rho_app =VES1D(R,thk,AB2);
%Inversion Parameter
nlayer = 3; 
nitr = 200; 
T = 5;
dec = 0.05;
%Generate random model
LBR = [1 1 1];
UBR = [100 500 50];
LBT = [1 1];
UBT = [25 50];
rho1(1 , :) = LBR + rand*(UBR - LBR);
thick1(1, :) = LBT + rand*(UBT - LBT);

%Calculate respon and misfit for each model
[rho_app_1] = VES1D(rho1(1,:),thick1(1,:),AB2);
app_rho1(1,:) = rho_app_1;
[misfit1] = misfit_VES(rho_app,app_rho1(1,:));
E1 = misfit1;

%Inversion process
for itr = 1 : nitr
    rho2(1 , :) = LBR + rand*(UBR - LBR);
    thick2(1, :) = LBT + rand*(UBT - LBT);
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
    T = T*(1-dec);
    Temperature(itr) = T;
end

%Data vizualization 
r_plot = [0, R];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho1];
Depth_mod = [0,cumsum(thick1),max(thick1)*100];
figure(1)
subplot(1,6,[1 3])
loglog(AB2,rho_app_new,'r',AB2,rho_app,'ob','MarkerSize',6,'LineWidth',2.5);
axis([1 10^3 1 10^3]);
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
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','Northeast');
axis([1 10^4 0 100]);
xlabel('Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
subtitle(['\rho_{1} = ',num2str(rho1(1)),' || \rho_{2} = ',num2str(rho1(2)),' || \rho_{3} = ',num2str(rho1(3)),' || thick_{1} = ',num2str(thick1(1)),' || thick_{2} = ',num2str(thick1(2))],'FontWeight','bold')
set(gca,'YDir','Reverse');
set(gca, 'XScale', 'log');
set(gcf, 'Position', get(0, 'Screensize'));
grid on


figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RSME','FontSize',10,'FontWeight','Bold');
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

