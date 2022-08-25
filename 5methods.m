%Minoo Ahmadi - 810897032
clc,clear all
close all
% Initialization
N=1e3;
A = data;
m = 10;
n = 20;

% 1: Equal Lamp Powers

%Graphically Using Plot
g_s = linspace(0,1,20)';
gamma = repmat(g_s, 1,10);
res = max(abs(log(data*transpose(gamma))));
plot(g_s, res, 'color',"#16a085",'linewidth',2.5)
title("Equal Lamp Powers")
grid on
print('equal_lamp.png','-dpng','-r300')
p_equal_lamp=min(res)*ones(m,1);

%Numerically
p = logspace(-3,0,N);
f = zeros(size(p));
for k=1:N
    f(k) = max(abs(log(A*p(k)*ones(m,1))));
end
[values_equal_lamp,imin] = min(f);




% 2: Least-squares with Rounding
p_least_squares_sat = A\ones(n,1);
p_least_squares_sat = max(p_least_squares_sat,0);
p_least_squares_sat = min(p_least_squares_sat,1);
p_least_squares_rnd = max(abs(log(A*p_least_squares_sat)));

% 3: Weighted Least-squares
rhos = linspace(1e-3,1,N);
crt = [];
for j=1:N
    p = [A; sqrt(rhos(j))*eye(m)]\[ones(n,1); sqrt(rhos(j))*0.5*ones(m,1)];
    crt = [ crt norm(p-0.5,inf) ];
end
indices = find(crt <= 0.5);
rho = rhos(indices(1)); 
p_least_squares_reg = [A; sqrt(rho)*eye(m)]\[ones(n,1); sqrt(rho)*0.5*ones(m,1)];
values_least_squares_reg = max(abs(log(A*p_least_squares_reg)));

% 4: Chebyshev Approximation
cvx_begin
    variable p_chbshv(m)
    minimize(norm(A*p_chbshv-1, inf))
    subject to
        p_chbshv >= 0
        p_chbshv <= 1
cvx_end
val_p_chbshv = max(abs(log(A*p_chbshv)));

% 5: Exact Solution Using CONVEX Lib:
cvx_begin
    variable p_exct(m)
    minimize(max([A*p_exct; inv_pos(A*p_exct)]))
    subject to
        p_exct >= 0
        p_exct <= 1
cvx_end

val_exct = max(abs(log(A*p_exct)));

% Answers to Problem Using 5 deifferent methods
[p_equal_lamp p_least_squares_sat p_least_squares_reg p_chbshv p_exct]
[values_equal_lamp p_least_squares_rnd values_least_squares_reg val_p_chbshv val_exct]
