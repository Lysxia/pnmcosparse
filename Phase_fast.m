function [] = Phase_fast(dir,d_,delta_,rho_,e_,n_,cores)

str2int = @(s) sscanf(s,'%d');

d = str2int(d_);
delta = str2int(delta_);
rho = str2int(rho_);
e = str2int(e_);
n = str2int(n_);
cores = str2int(cores);

if (cores > 1 && matlabpool('size') == 0)
    matlabpool('open',cores);
end

if ~isdir(dir)
    disp('Not an existing directory');
    return
end

%% Phase Diagram Drawing

param_str = strcat('dim',d_,'d',delta_,'r',rho_,'e',e_,'noise',n_);
Omega = fdamatrix(d);

%% ACoSaMP a=1 lambda=1

disp({'ACoSaMP'});

a = 1;
lambda = 1;

param = @(y,M,l) params_analysis(d,y,M,Omega,l,lambda,a,'projector');

[D_Facosamp,var_Facosamp] = csa_phasediag( d,delta,rho,@asp,param,n,e);

name_acosamp = strcat('Facosamp','_',param_str);
save(strcat(dir,'/',name_acosamp,'.mat'),'D_Facosamp','var_Facosamp');

%% ASP a=1 lambda=1

disp({'ASP'});

a = 1;
lambda = 1;

param = @(y,M,l) params_analysis(d,y,M,Omega,l,lambda,a,'estimate');

[D_Fasp,var_Fasp] = csa_phasediag( d,delta,rho,@asp,param,n,e);

name_asp = strcat('Fasp','_',param_str);
save(strcat(dir,'/',name_asp,'.mat'),'D_Fasp','var_Fasp');

%% AIHT (adaptive gradient step)

disp({'AIHT'});


param = @(y,M,l) params_analysis(d,y,M,Omega,l,lambda,a,'projector');

[D_Faiht,var_Faiht] = csa_phasediag( d,delta,rho,@aiht,param,n,e);

name_aiht = strcat('Faiht','_',param_str);
save(strcat(dir,'/',name_aiht,'.mat'),'D_Faiht','var_Faiht');

%% AHTP (adaptive gradient step) lambda = 1

disp({'AHTP'});

lambda = 1;

param = @(y,M,l) params_analysis(d,y,M,Omega,l,lambda,a,'estimate');

[D_Fahtp,var_Fahtp] = csa_phasediag( d,delta,rho,@aiht,param,n,e);

name_ahtp = strcat('Fahtp','_',param_str);
save(strcat(dir,'/',name_ahtp,'.mat'),'D_Fahtp','var_Fahtp');

if cores > 1
    matlabpool('close');
end