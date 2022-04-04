function [misfit] = misfit_VES(rho_obs, rho_cal)
    ls = length(rho_obs);
    for j = 1 : ls
        m(j) = ((rho_cal(j) - rho_obs(j)))^2;
    end
    misfit = sqrt((1/ls)*sum(m));
end
