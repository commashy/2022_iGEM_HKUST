%% This ODE represents the HIV model in Section 4.2
function dzdt=ODEmodel(t,y,X,run_num)

%% PARAMETERS %%
Parameter_settings_EFAST;

k1=X(run_num,1);
k2=X(run_num,2);
k3=X(run_num,3);
kd3=X(run_num,4);
k4=X(run_num,5);
k5=X(run_num,6);
k6=X(run_num,7);
kd6=X(run_num,8);
k7=X(run_num,9);
k8=X(run_num,10);
gam1=X(run_num,11);
gam2=X(run_num,12);
n=X(run_num,13);
m=X(run_num,14);
d1=X(run_num,15);
d2=X(run_num,16);
d3=X(run_num,17);
d4=X(run_num,18);
d5=X(run_num,19);
d6=X(run_num,20);
d7=X(run_num,21);
d8=X(run_num,22);
dummy=X(run_num,9);

% INITIAL CONDITION FOR THE ODE MODEL
%y(1) = dia = 10000;
%h2o2=0;
%oxyrn=0;
%oxyra=0;
%poxy=0.0001;
%poxy_oxyra=0;
%mgfp=0;
%gfp=0;
%pkat=0.0001;
%pkat_oxyra=0;
%mrfp=0;
%rfp=0;

% ODEs
ddiadt = 0.01 -k1*y(1) - d5*y(1);
dh2o2dt = k1*y(1) -k2*y(3)*y(2) + gam2 - d6*y(2);
doxyrndt = gam1 - k2*y(3)*y(2) - d7*y(3);
doxyradt = k2*y(3)*y(2) - d8*y(4);
dpoxydt = - k3*y(5)*y(4);
dpoxy_oxyradt = k3*y(5)*y(4); 
dmgfpdt = n - k4*y(7) - d1*y(7);
dgfpdt = k5*y(7) - d2*y(8);
dpkatdt = - k6*y(9)*y(4)^m;
dpkat_oxyradt = k6*y(9)*y(4)^m;
dmrfpdt = m - k7*y(11)  - d3*y(11);
drfpdt = k8*y(11) - d4*y(12);
        
dzdt = [ddiadt; dh2o2dt; doxyrndt; doxyradt;
        dpoxydt; dpoxy_oxyradt; dmgfpdt; dgfpdt;
        dpkatdt; dpkat_oxyradt; dmrfpdt; drfpdt];


% [T] CD4+ uninfected: Tsource + Tprolif - Tinf
%Tsource = s - muT*y(1);
%Tprolif = r*y(1)*(1-(y(1)+y(2)+y(3))/Tmax);
%Tinf = k1*y(1)*y(4);

% [T1] CD4+ latently infected: Tinf - T1death - T1inf
%T1death = muT*y(2);
%T1inf = k2*y(2);

% [T2] CD4+ actively infected: T1inf - T2death
%T2death = mub*y(3);

% [V] Free infectious virus: Vrelease - Tinf - Vdeath
%Vrelease = N*T2death;
%Vdeath = muV*y(4);

%dydt = [Tsource + Tprolif - Tinf;
%        Tinf - T1death - T1inf;
%        T1inf - T2death;
%        Vrelease - Tinf - Vdeath];