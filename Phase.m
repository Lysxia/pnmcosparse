function [] = Phase(dir,d_,delta_,rho_,e_,n_,cores)

str2int = @(s) sscanf(s,'%d');

d = str2int(d_);
delta = str2int(delta_);
rho = str2int(rho_);
e = str2int(e_);
n = str2int(n_);
cores = str2int(cores);

matlabpool('open',cores);

if ~isdir(dir)
    disp('Not an existing directory');
    return
end

%% Phase Diagram Drawing

param_str = strcat('dim',d_,'d',delta_,'r',rho_,'e',e_,'noise',n_);

%% ACoSaMP a=1 lambda=1

disp({'ACoSaMP'});

a = 1;
lambda = 1;

param = @(y,M,l) params_1DFD(d,y,M,l,lambda,a,'projector');

[D_acosamp,var_acosamp] = csa_phasediag( d,delta,rho,@asp,param,n,e);

name_acosamp = strcat('acosamp','_',param_str);
save(strcat(dir,'/',name_acosamp,'.mat'),'D_acosamp','var_acosamp');

%% ASP a=1 lambda=1

disp({'ASP'});

a = 1;
lambda = 1;

param = @(y,M,l) params_1DFD(d,y,M,l,lambda,a,'estimate');

[D_asp,var_asp] = csa_phasediag( d,delta,rho,@asp,param,n,e);

name_asp = strcat('asp','_',param_str);
save(strcat(dir,'/',name_asp,'.mat'),'D_asp','var_asp');

%% AIHT (adaptive gradient step)

disp({'AIHT'});


param = @(y,M,l) params_1DFD(d,y,M,l,lambda,a,'projector');

[D_aiht,var_aiht] = csa_phasediag( d,delta,rho,@aiht,param,n,e);

name_aiht = strcat('aiht','_',param_str);
save(strcat(dir,'/',name_aiht,'.mat'),'D_aiht','var_aiht');

%% AHTP (adaptive gradient step) lambda = 1

disp({'AHTP'});

lambda = 1;

param = @(y,M,l) params_1DFD(d,y,M,l,lambda,a,'estimate');

[D_ahtp,var_ahtp] = csa_phasediag( d,delta,rho,@aiht,param,n,e);

name_ahtp = strcat('ahtp','_',param_str);
save(strcat(dir,'/',name_ahtp,'.mat'),'D_ahtp','var_ahtp');

matlabpool('close');