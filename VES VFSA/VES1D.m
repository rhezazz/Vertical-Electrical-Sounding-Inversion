%Program pemodelan kedepan kurva sounding resistivitas 1-D dengan digital filter
% Reference I: Ekinci Y L and Demirci A. 2008. J. of Appl. Sciences 8 (22): 4070-4078
function rho_app = VES1D(res,thk,AB2)
indmin = -2;
indmax = 10;
ind = indmin:indmax;
psi = 4.438;
ls = length(AB2);
koef_filter = [225, 336, 2525, 8183, 16448, -27841, 13396, -4390,1605, -746, 416, -262, 105]/10000;
for i = 1 : ls
    lambda = (1/AB2(i))*exp(ind*log(10)/psi);
    lay = length(res);
    for j = 1 : length(ind)
        T = res(lay);
        w = lay;
        while w>1
            w = w-1;
            aa = tanh(thk(w)*lambda(j));
            T = (T+res(w)*aa)/(1+T*aa/res(w));
        end
        Ti(j) = T;
    end
    rho_app(i) = koef_filter*Ti';
end
end