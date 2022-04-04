%Flower Pollination Algorithm for VES Data Inversion
%Mohammad Rheza Zamani
%Reference : Xin-She Yang. 2014. Nature-Inspired Optimization Algorithms (1st. ed.). Elsevier Science Publishers B. V., NLD.
tic;
clear all;
clc;
%Synthetic Model
k = 0:17;
AB2 = 10.^(k/6);
res = [500 250 1000 500];
thk = [5 25 50];
rho_app = VES1D(res,thk,AB2);
%Inversion Parameter
npop = 50; 
nlayer = 4; 
nitr = 5000; 
prp = 0.8;
%Generate random model
LBR = [1 1 1 1];
UBR = [2000 2000 2000 2000];
LBT = [1 1 1];
UBT = [50 50 50];
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
%Find global best
idx = find(E ==min(E));
rho_best = rho(idx(1),:);
thick_best = thick(idx(1),:);
%Inversion process
for iitr =  1 : nitr
    for i = 1 : npop
        %Global Pollination
        if rand < prp
            n = 1;
            m1 = nlayer;
            beta = 1.5;
            [Lr]=levy(n,m1,beta);
            for n = 1 : nlayer
                npopr(1,n)=rho(i,n)+(0.1.*Lr(n).*(rho(i,n)-rho_best(n)));
                if npopr(1,n) < LBR(n);
                     npopr(1,n) = LBR(n);
                end
                if npopr(1,n) > UBR(n);
                     npopr(1,n) = UBR(n);
                end
            end
             m2 = nlayer -1;
            [Lt]=levy(n,m2,beta);
            for n = 1 : (nlayer-1)
                npopt(1,n)=thick(i,n)+(0.1.*Lt(n).*(thick(i,n)-thick_best(n)));
                if npopt(1,n) < LBT(n);
                     npopt(1,n) = LBT(n) ;
                end
                if npopt(1,n) > UBT(n);
                    npopt(1,n) = UBT(n);
                 end
            end
        else
            %Local Pollination
            epsilon=rand;
            JK=randperm(npop);
            for n = 1 : nlayer
                npopr(1,n)=rho(i,n)+epsilon*(rho(JK(1),n)-rho(JK(2),n));
                 if npopr(1,n) < LBR(n);
                     npopr(1,n) = LBR(n);
                 end
                 if npopr(1,n) > UBR(n);
                    npopr(1,n) = UBR(n);
                 end
            end
            for n = 1 : (nlayer-1)
                npopt(1,n)=thick(i,n)+epsilon*(thick(JK(1),n)-thick(JK(2),n));
                 if npopt(1,n) < LBT(n);
                     npopt(1,n) = LBT(n);
                 end
                 if npopt(1,n) > UBT(n);
                     npopt(1,n) = UBT(n);
                end
            end
        end
        %Calculate respon and misfit for each new model
        rho_new=npopr;
        thick_new=npopt;
        [rho_app_new]= VES1D(rho_new,thick_new,AB2);
        [err] = misfit_VES(rho_app,rho_app_new);
        %Update model for each model
        if err<E(i)
            rho(i,:) = rho_new(1,:);
            thick(i,:) = thick_new(1,:);
            app_res(i,:) = rho_app_new(1,:);
            E(i) = err;
        end
    end
    %Update model for each iteration
    Emin = 100;
        for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            rho_best = rho(ipop,:);
            thick_best = thick(ipop,:);
            app_model = app_res(ipop,:);
        end
    end
    Egen(iitr)=Emin;
end

%Data vizualization 
r_plot = [0, res];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,rho_best];
Depth_mod = [0,cumsum(thick_best),max(thick_best)*100];

figure(1)
subplot(1,6,[1 3])
loglog(AB2,app_model,'r',AB2,rho_app,'ob','MarkerSize',6,'LineWidth',2.5);
axis([1 10^3 1 10^4]);
legend({'Observed Data','Calculated Data'},'Color','none','FontWeight','Bold');
xlabel('AB/2 (m)','FontSize',8,'FontWeight','Bold');
ylabel('App. Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}Respon  || Misfit : ', num2str(Egen(iitr)),' || iteration : ', num2str(iitr)]);
grid on
subplot(1,6,[5 6])
stairs(r_plot,t_plot,'--r','Linewidth',2.5);
hold on
stairs(r_mod,Depth_mod,'-b','Linewidth',1.5);
hold off
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','Southwest');
axis([1 10^4 0 200]);
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
ylabel('RSME','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
grid on
time = toc;


%Levy Function
function [z] = levy(n,m,beta)
    num = gamma(1+beta)*sin(pi*beta/2);
    
    den = gamma((1+beta)/2)*beta*2^((beta-1)/2);

    sigma_u = (num/den)^(1/beta);

    u = normrnd(0,sigma_u^2,n,m); 
    
    v = normrnd(0,1,n,m);

    z = u./(abs(v).^(1/beta));
end