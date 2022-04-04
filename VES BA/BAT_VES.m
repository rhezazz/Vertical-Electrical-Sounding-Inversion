%Bat algorithm for VES Data Inversion
%Mohammad Rheza Zamani
%Reference : Xin-She Yang. 2014. Nature-Inspired Optimization Algorithms (1st. ed.). Elsevier Science Publishers B. V., NLD.
clc;
clear all;
tic
%Synthetic Model
n = 0:17;
AB2 = 10.^(n/6);
res = [500 1000 250 1500];
thk = [4 15 25];
rho_app = VES1D(res,thk,AB2);
nlayer = length(res);
alpha = 0.9;
gamma = 0.9;
%Inversion Parameter
npop = 200;
niter = 1000;
%Minimum and maximum frequency 
fmin = 0;
fmax = 2;
%Velocity Initialization
%Calculate initial velocity
for i = 1 : npop
    for imod =  1 : nlayer
        v_rho(i,imod) = 0;
    end
    for imod = 1 : nlayer -1
        v_thk(i,imod) = 0;
    end
end
%Minimum and maximum loudness
Amin = 1;
Amax = 2;
%Loudness Initialization
%Calculate initial loudness
for i = 1 : npop
    for imod = 1 : nlayer
        A_rho(i,imod) = Amin + rand*(Amax - Amin);
    end
    for imod =  1 : nlayer-1
        A_thk(i,imod) = Amin + rand*(Amax - Amin);
    end
end
%Minimum and maximum pulsa rate
rmin = 0;
rmax = 1;
%Pulsa rate Initialization
%Calculate initial pulse rate
for i = 1 : npop
    for imod = 1 :  nlayer
        r_rho(i,imod) = rmin + rand*(rmax-rmin);
    end
    for  imod = 1 : nlayer-1
        r_thk(i,imod) = rmin + rand*(rmax-rmin);
    end
end
%Generate random model
rhomin = [1 1 1 1];
rhomax = [1500 1500 1500 1500];
thkmin = [1 1 1];
thkmax = [30 30 30];
for i = 1 : npop
    rho(i , :) = rhomin + rand*(rhomax - rhomin);
    thick(i, :) = thkmin + rand*(thkmax - thkmin);
end
%Calculate respon and misfit for each model
for i = 1 : npop
    [res_semu] = VES1D(rho(i,:),thick(i,:),AB2);
    app_res(i,:) = res_semu;
    
    [misfit] = misfit_VES(rho_app,app_res(i,:));
    E(i) = misfit;
end

%Find global best
idx = find(E ==min(E));
G_best_rho = rho(idx(1),:);
G_best_thick = thick(idx(1),:);
%Inversion process
for itr =  1 : niter
    for i = 1 :  npop
        %Frequency initialization
        for imod = 1 : nlayer
            f_rho(imod) = fmin + rand*(fmin-fmax);
        end
        for imod = 1 : nlayer-1
            f_thk(imod) = fmin + rand*(fmin-fmax);
        end
        %Update velocity and position
        for imod  = 1 : nlayer
            v_rho(1,imod) = v_rho(i,imod) + (rho(i,imod)-G_best_rho(imod))*f_rho(imod);
            rho_baru(1,imod) = rho(i,imod) + v_rho(1,imod);
            if rho_baru(1,imod) > rhomax(imod)
                rho_baru(1,imod) = rhomax(imod);
            end
            if rho_baru(1,imod) < rhomin(imod)
                rho_baru(1,imod) = rhomin(imod);
            end
        end
       for imod  = 1 : nlayer-1
            v_thk(1,imod) = v_thk(i,imod) + (thick(i,imod)-G_best_thick(imod))*f_thk(imod);
            thick_baru(1,imod) = thick(i,imod) +v_thk(1,imod);
            if thick_baru(1,imod) > thkmax(imod)
                thick_baru(1,imod) = thkmax(imod);
            end
            if thick_baru(1,imod) < thkmin(imod)
                thick_baru(1,imod) = thkmin(imod);
            end
       end
       random = rand;
       if random > r_rho(i) 
          for imod = 1 : nlayer
              rho_baru(1,imod) =G_best_rho(imod)+mean(A_rho(i,:))*(-1+2*rand);
          end
           if rho_baru(1,imod) > rhomax(imod)
                rho_baru(1,imod) = rhomax(imod);
            end
            if rho_baru(1,imod) < rhomin(imod)
                rho_baru(1,imod) = rhomin(imod);
            end
       end
       if random > r_thk(i)
           for imod = 1 : nlayer-1
               thick_baru(1,imod) = G_best_thick(imod) + mean(A_thk(i,:))*(-1+2*rand);
           end
             if thick_baru(1,imod) > thkmax(imod)
                thick_baru(1,imod) = thkmax(imod);
            end
            if thick_baru(1,imod) < thkmin(imod)
                thick_baru(1,imod) = thkmin(imod);
            end
       end
       %Calculate respon and misfit for each new model
       [app_res_baru] = VES1D(rho_baru,thick_baru,AB2);
       [E_baru] = misfit_VES(rho_app,app_res_baru);
       %Update model for each model
       if E_baru < E(i) && rand < A_rho(i) && rand < A_thk(i) 
            rho(i,:) = rho_baru(1,:);
            thick(i,:) = thick_baru(1,:);
            E(i) = E_baru;
            app_res(i,:) = app_res_baru(1,:);
       end
    end
    %Update model for each iteration
    Emin = 1000;
    for i = 1 : npop
      if E(i)< Emin
         Emin = E(i);
         G_best_rho  = rho(i,:);
         G_best_thick  = thick(i,:);
         rho_app_mod = app_res(i,:);
      end
    end
     A_rho(i) = A_rho(i)*alpha;
     A_thk(i) = A_thk(i)*alpha;
     r_rho(i) = r_rho(i)*(1-exp(gamma*itr));
     r_thk(i) = r_thk(i)*(1-exp(gamma*itr));
    Egen(itr)=Emin;
end
toc

%Data vizualization 
r_plot = [0, res];
t_plot = [0,cumsum(thk),max(thk)*100];
r_mod = [0,G_best_rho];
Depth_mod = [0,cumsum(G_best_thick),max(G_best_thick)*100];

figure(1)
subplot(2,2,[1 3])
loglog(AB2,rho_app_mod,'r',AB2,rho_app,'ob','MarkerSize',5,'LineWidth',3);
axis([1 10^3 1 10^4]);
legend({'Calculated Data','Observed Data'},'Color','none','FontWeight','Bold');
xlabel('AB/2 (m)','FontSize',8,'FontWeight','Bold');
ylabel('App. Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
title(['\bf \fontsize{10}\fontname{Times}Respon  || Misfit : ', num2str(Egen(itr)),' || iteration : ', num2str(itr)]);
grid on
subplot(2,2,[2 4])
stairs(r_plot,t_plot,'--r','Linewidth',3);
hold on
stairs(r_mod,Depth_mod,'-b','Linewidth',2);
hold off
legend({'Synthetic Model','Inversion Model'},'Color','none','FontWeight','Bold','Location','southwest');
axis([1 10^4 0 100]);
xlabel('Resistivity (Ohm.m)','FontSize',8,'FontWeight','Bold');
ylabel('Depth (m)','FontSize',8,'FontWeight','Bold');
title('\bf \fontsize{10} Model');
set(gca,'YDir','Reverse');
set(gca, 'XScale', 'log');
set(gcf, 'Position', get(0, 'Screensize'));
grid on

figure(2)
plot(1:niter,Egen,'r','Linewidth',2)
xlabel('Iteration Number','FontSize',10,'FontWeight','Bold');
ylabel('RMSE','FontSize',10,'FontWeight','Bold');
set(gcf, 'Position', get(0, 'Screensize'));
grid on