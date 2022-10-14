function Fisherly_Modeling_Zero_new(tf, param)

%Fisherly_Modeling_Zero_new(10000, [0.01,0,0,0,1,1,0,0,0,0,0])

% Initial Concentrations
a0=param(1); %Diamine
s0=param(2); %H2O2
d0=param(3); %OxyRn
f0=param(4); %OxyRa
g0=param(5); %Pkat
h0=param(6); %PoxyS
j0=param(7); %GFP
k0=param(8); %RFP
z0=param(9); %mRNAGFP
x0=param(10); %mRNARFP
N0 = [a0;s0;d0;f0;g0;h0;j0;k0;z0;x0];


% Rate Constants
k1 = 0.00009627; 
k2 = 0.00042; 
k4 = 0.00001;
k5 = 0.1;
k6 = 0.3;
alpha = 0.5; 
beta = 0.1;
k = 0.0001; 
L = 0.001; 
k3 = 0.01; 


m1 = 0;
m2 = 0.0001;

d2 = 0.00000802;
d3 = 0.001;
d4 = 0.001;
d7 = 0.0021; 
d8 = 0.000138; 
d5 = 0.1; 
d6 = 0.1; 
d1 = 0.01; 


%Running the Built-in ODE Solver
[tt,NN] = ode45(@f,[0,tf],N0);

AA = NN(:,1);
SS = NN(:,2);
DD = NN(:,3);
FF = NN(:,4);
GG = NN(:,5);
HH = NN(:,6);
JJ = NN(:,7);
KK = NN(:,8);
ZZ = NN(:,9); 
XX = NN(:,10);


%TARGET GRAPH
figure()
plot(tt, JJ, 'green', tt, KK, 'red');
title('Target Graph');
%xline(6000, '-','Point of Inspection');
xlabel('Time (Sec)');
ylabel('Fluorescence Concentration (Molar)');
xline(6000, '-','Point of Inspection')
ylim([0 7]);
xlim([0 inf]);

    function dYdt = f(t,Y)
        Dia=Y(1);
        H2O2=Y(2);
        OxyRn=Y(3);
        OxyRa=Y(4);
        Pkat=Y(5);
        PoxyS = Y(6);
        GFP =Y(7);
        RFP=Y(8);
        mRNAGFP =Y(9);
        mRNARFP=Y(10);  

        %ODEs that govern the concentrations
        dDiadt = m1 - k1.*Dia - d1.*Dia;
        dH2O2dt = k1.*Dia - k2.*OxyRn.*H2O2 - d2.*H2O2;
        dOxyRndt = m2 - k2.*OxyRn.*H2O2 - d3.*OxyRn;
        dOxyRadt = k2.*OxyRn.*H2O2 - d4.*OxyRa;
        dPoxySdt = - k3 .* OxyRa .* PoxyS; 
	    dPkatdt = - k4.*OxyRa.*Pkat;
        dmRNAGFPdt = alpha.*OxyRa./(k+OxyRa) - d5.*mRNAGFP; 
        dmRNARFPdt = beta*OxyRa./(L+OxyRa) - d6.*mRNARFP;

	    dGFPdt = k5.*mRNAGFP - d7.*GFP;
	    dRFPdt = k6.*mRNARFP - d8.*RFP;

        dYdt=[dDiadt;dH2O2dt;dOxyRndt;dOxyRadt;dPkatdt;dPoxySdt;dGFPdt;dRFPdt;dmRNAGFPdt;dmRNARFPdt];



    end


end