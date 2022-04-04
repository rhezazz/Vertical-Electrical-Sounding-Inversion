%Fireflly Algorithm for VES Data Inversion
%Mohammad Rheza Zamani
%Reference : Xin-She Yang. 2014. Nature-Inspired Optimization Algorithms (1st. ed.). Elsevier Science Publishers B. V., NLD.
tic
clear all;
clc;
%Synthetic Model
n = 0:17;
AB2 = 10.^(n/6);
res = [100 20 50];
thk = [5 10];
rho_app = VES1D(res,thk,AB2);
%Inversion Parameter
npop = 100; 
nlayer = 3; 
nitr = 500; 
%Generate random model
LBR = [1 1 1];
UBR = [500 100 250];
LBT = [1 1];
UBT = [25 50];
alpha = 0.2;
betha0 = 1;
gamma = 0.8;
damp = 0.99;
for ipop = 1 : npop
    rho(ipop , :) = LBR + rand*(UBR - LBR);
    thick(ipop, :) = LBT + rand*(UBT - LBT);
end
%Calculate respon and misfit for each model
for ipop = 1 : npop 
    [res_semu] = VES1D(rho(ipop,:),thick(ipop,:),AB2);
    app_res(ipop,:) = res_semu;
    
    [misfit] = misfit_VES(rho_app,app_res(ipop,:));
    E(ipop) = misfit;
end

%Inversion process
for itr = 1 : nitr
    for i =  1 : npop
        j = randi(npop,1);
        while i == j
            j = randi(npop,1);
        end
        if E(i)<E(j)
            %Calculated Distance
            dr = norm((rho(i,:)-rho(j,:)));
            dt = norm((thick(i,:)-thick(j,:)));
            %Calculated new model with determine a new position
            %Random vector position for resistivity model
            for n1 = 1 : nlayer
                
                rho_baru(1,n1) = rho(i,n1) +betha0.*exp(-gamma*(dr)^2).*(rho(j,n1)-rho(i,n1))+ (alpha*(rand-0.5)*abs((UBR(n1)-LBR(n1))));
                if rho_baru(1,n1) < LBR(n1);
                     rho_baru(1,n1) = LBR(n1);
                end
                if rho_baru(1,n1) > UBR(n1);
                     rho_baru(1,n1) = UBR(n1);
                end
            end
            %Random vector position for thick model
            for n2 = 1 : (nlayer-1)
                thk_baru(1,n2) = thick(i,n2) +betha0.*exp(-gamma*(dt)^2).*(thick(j,n2)-thick(i,n2))+ (alpha*(rand-0.5)*abs((UBT(n2)-LBT(n2))));
                if thk_baru(1,n2) < LBT(n2);
                     thk_baru(1,n2) = LBT(n2);
                end
                if thk_baru(1,n2) > UBT(n2);
                    thk_baru(1,n2) = UBT(n2);
                end
            end 
        else
            rho_baru(1,:) = rho(i,:);
            thk_baru(1,:) = thick(i,:);
        end
        %Calculate respon and misfit for each new model
        [res_app_baru]= VES1D(rho_baru,thk_baru,AB2);
        [err] = misfit_VES(rho_app,res_app_baru);
         %Update model for each model
        if err<E(i)
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thk_baru(1,:);
            app_res(i,:) = res_app_baru(1,:);
            E(i) = err;
        end
      
    end
    %Update model for parasitism phase
     Emin = 100;
     for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            rho_model = rho(ipop,:);
            thk_model = thick(ipop,:);
            app_model = app_res(ipop,:);
        end
    end
    Egen(itr)=Emin;
    alpha = alpha*damp;
end
time = toc;

%Data vizualization 
r_plot = [0, res];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho_model];
Depth_mod = [0,cumsum(thk_model),max(thk_model)*100];

figure(1)
subplot(1,6,[1 3])
loglog(AB2,app_model,'r',AB2,rho_app,'ob','MarkerSize',6,'LineWidth',2.5);
axis([1 10^3 1 10^3]);
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('AB/2 (m)','FontSize',8,'FontWeight','Bold');
ylabel('App. Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}Respon  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
grid on
subplot(1,6,[5 6])
stairs(r_plot,t_plot,'--r','Linewidth',2);
hold on
stairs(r_mod,Depth_mod,'-b','Linewidth',2);
hold off
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','southwest');
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
grid on
