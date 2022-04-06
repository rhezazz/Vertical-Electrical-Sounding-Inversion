%Program pemodelan inversi kurva sounding resistivitas 1-D dengan
%menggunakan algoritma PSO (Particle Swarm Optimization)
%Mohammad Rheza Zamani
clear all;
clc;
%Input data sintetik
n = 0:17;
AB2 = 10.^(n/6);
res = [50 100 20];
thk = [10 20];
%Pehitungan data sintetik
rho_app = VES1D(res,thk,AB2);
%Definisi ruang model 
npop = 100; %Jumlah dari model  
nlayer = 3; %Jumlah lapisan 
nitr = 1000; %Jumlah iterasi
%Paramter inversi
wmax = 0.9;
wmin = 0.5;
c1=2.05;
c2 =2.05;
%Batas atas dan bawah (Diatur 5 kali dari nilai parameter model)
%Batas bawah pencarian nilai resistivitas
rhomin = [1 1 1];
%Batas atas pencarian nilai resistivitas
rhomax = [1000 1000 1000];
%Batas bawah pencarian nilai ketebalan
thkmin = [1 1];
%Batas atas pencarian nilai resistivitas
thkmax = [50 50];
%Membuat model awal acak
for ipop = 1 : npop
    rho(ipop,:) = rhomin + rand*(rhomax - rhomin);
    thick(ipop,:) = thkmin + rand*(thkmax - thkmin);
end
%Hitung velocity awal dan posisi awal
for ipop = 1 : npop
    for imod =  1 : nlayer
        %v_rho(ipop,imod) = 0.5.*(min(rho(ipop,:)) + rand*(max(rho(ipop,:)) - min(rho(ipop,:))));
        v_rho(ipop,imod) = 0;
    end
    for imod = 1 : nlayer -1
        %v_thk(ipop,imod) = 0.5.*(min(thick(ipop,:))) + rand*(max(thick(ipop,:)) - min(thick(ipop,:)));
        v_thk(ipop,imod) = 0;
    end
end
%Menhitung masing - masing resistivitas model semu dan misfit dari model
for ipop = 1 : npop 
    [res_semu] = VES1D(rho(ipop,:),thick(ipop,:),AB2);
    app_res(ipop,:) = res_semu;
    
    [misfit] = misfit_VES(rho_app,app_res(ipop,:));
    E(ipop) = misfit;
end
%Global best
idx = find(E ==min(E));
G_best_rho = rho(idx(1),:);
G_best_thick = thick(idx(1),:);
%Inversi
for itr = 1 : nitr
    w = wmin+((wmax-wmin)/nitr)*itr;
    for i = 1 : npop
        %personal best
        P_best_rho = rho;
        P_best_thick = thick;
        %Membuat komponen kecepatan
        %Rho
        for n1 = 1 : nlayer
            v_rho(1,n1) = w.*v_rho(i,n1) + c1.*rand.*(P_best_rho(i,n1) - rho(i,n1))+ c2.*rand.*(G_best_rho(n1) - rho(i,n1));
            rho_baru(1,n1) = rho(i,n1)+ v_rho(1,n1);
        if rho_baru(1,n1)<rhomin(n1)
            rho_baru(1,n1) = rhomin(n1);
        end
        if rho_baru(1,n1)>rhomax(n1)
            rho_baru(1,n1) = rhomax(n1);
        end
        end
        %Ketebalan
        for n2 = 1 : (nlayer-1)
            v_thk(1,n2) = w.*v_thk(i,n2) + c1.*rand.*(P_best_thick(i,n2) - thick(i,n2))+ c2.*rand.*(G_best_thick(n2) - thick(i,n2));
            thk_baru(1,n2) = thick(i,n2)+ v_thk(1,n2);
          if thk_baru(1,n2)<thkmin(n2)
              thk_baru(1,n2) = thkmin(n2);
          end
          if thk_baru(1,n2)>thkmax(n2)
              thk_baru(1,n2) = thkmax(n2);
          end
        end
        %Update pesonal best
        [res_app_baru]= VES1D(rho_baru,thk_baru,AB2);
        [E_baru] = misfit_VES(rho_app,res_app_baru);
        if E_baru<E(i)
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thk_baru(1,:);
            app_res(i,:) = res_app_baru(1,:);
            E(i) = E_baru;
        end
    end
    Emin = 100;
    for ipop = 1 : npop
        if E(ipop)< Emin
            Emin = E(ipop);
            G_best_rho  = rho(ipop,:);
            G_best_thick  = thick(ipop,:);
            rho_app_mod = app_res(ipop,:);
        end
    end
    Egen(itr)=Emin;
end
%Ploting model
r_plot = [0, res];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,G_best_rho];
Depth_mod = [0,cumsum(G_best_thick),max(G_best_thick)*100];

figure(1)
subplot(1,6,[1 3])
loglog(AB2,rho_app_mod,'r',AB2,rho_app,'ob','MarkerSize',6,'LineWidth',2.5);
axis([1 10^3 1 10^3]);
legend({'Observed Data','Calculated Data'},'Color','none','FontWeight','Bold');
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
axis([1 10^4 0 100]);
xlabel('Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
% subtitle(['\rho_{1} = ',num2str(G_best_rho(1)),' || \rho_{2} = ',num2str(G_best_rho(2)),' || \rho_{3} = ',num2str(G_best_rho(3)),' || thick_{1} = ',num2str(G_best_thick(1)),' || thick_{2} = ',num2str(G_best_thick(2))],'FontWeight','bold')
set(gca,'YDir','Reverse');
set(gca, 'XScale', 'log');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

%plot misfit
figure(2)
plot(1:nitr,Egen,'r','Linewidth',1.5)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RMSE','FontSize',10,'FontWeight','Bold');
title('\bf \fontsize{12} Grafik Misfit ');
set(gcf, 'Position', get(0, 'Screensize'));
grid on