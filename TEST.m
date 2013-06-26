%% TEST

d = 200;
m = 200;
% p = d-1
l = d-150;
l_target = l;

a = 0.8;
lambda = 1;

iter = 50;
threshold = 0.001;

x = rand_1DFD_cs(d,d-1-l);
M = rand_sampling_matrix(m,d);
y = M*x+0.1*(2*rand(m,1)-1);

x0 = zeros(d,1);
%x1 = x+0.5*rand(d,1);

D = norm(y-M*x0);

if strcmp(strat,'acosamp') || strcmp(strat,'aiht')
    projector = 'projector';
else
    projector = 'estimate';
end

params_FD = params_1DFD(d,y,M,l_target,a,lambda,projector);
params_generic = params_analysis(d,y,M,fdamatrix(d),l_target,a,lambda,projector);

%% Test 1DFD
Solve;