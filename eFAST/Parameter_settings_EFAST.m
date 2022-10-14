%% PARAMETER INITIALIZATION
%% the range of parameters is set to be +-20%
pmin=[
    0.000077 
0.000336 
0.008000 
0.003200 
0.008000 
0.080000 
0.000008 
0.005120 
0.080000 
0.240000 
0.000080 
0.000640 
0.400000 
0.010400 
0.080000 
0.001680 
0.008000 
0.000110 
0.008000 
0.000006 
0.000800 
0.000800 
1]; % dummy

pmax=[0.000116 
0.000504 
0.012000 
0.004800 
0.012000 
0.120000 
0.000012 
0.007680 
0.120000 
0.360000 
0.000120 
0.000960 
0.600000 
0.015600 
0.120000 
0.002520 
0.012000 
0.000166 
0.012000 
0.000010 
0.001200 
0.001200 
2]; % dummy

% Parameter Labels 
efast_var={'k_1','k_2','k_3' ...
    'kd_3','k_4','k_5','k_6','kd_6','k_7','k_8','gam1','gam2','n','m', 'd_1','d_2','d_3','d_4','d_5','d_6','d_7','d_8','dummy'}';%,

% PARAMETER BASELINE VALUES
 k1	=	0.00009627	;
k2	=	0.00042	;
k3	=	0.01	;
kd3	=	0.004	;
k4	=	0.01	;
k5	=	0.1	;
k6	=	0.00001	;
kd6	=	0.0064	;
k7	=	0.1	;
k8	=	0.3	;
gam1	=	0.0001	;
gam2	=	0.0008	;
n	=	0.5	;
m	=	0.013	;
d1	=	0.1	;
d2	=	0.0021	;
d3	=	0.01	;
d4	=	0.000138	;
d5	=	0.01	;
d6	=	0.00000802	;
d7	=	0.001	;
d8	=	0.001	;
 dummy = 1.5;
%% TIME SPAN OF THE SIMULATION
t_end=4000; % length of the simulations
tspan=(0:1:t_end);   % time points where the output is calculated
time_points=[2000 4000]; % time points of interest for the US analysis

% INITIAL CONDITION FOR THE ODE MODEL
dia=0.01;
h2o2=0;
oxyrn=0;
oxyra=0;
poxy=0.0001;
poxy_oxyra=0;
mgfp=0;
gfp=0;
pkat=0.0001;
pkat_oxyra=0;
mrfp=0;
rfp=0;

%T0=1e3;
%T1=0;
%T2=0;
%V=1e-3;

y0=[dia, h2o2, oxyrn, oxyra, poxy, poxy_oxyra, mgfp, gfp, pkat, pkat_oxyra, mrfp, rfp];
%y0=[T0,T1,T2,V];
%y0=[T0,T1,T2,V];

% Variables Labels
y_var_label={'Dia', 'H2O2', 'OxyRn', 'OxyRa', 'Pkat', 'PoxyS', 'GFP', 'RFP', 'mRNAGFP', 'mRNARFP'};
%y_var_label={'dia','T*','T**','V'};
