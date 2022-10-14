function [GFPt, RFPd] = Fisherly_Modeling(tf, param)

% Initial Concentrations
a0=param(1); %Diamine
s0=param(2); %H2O2
d0=param(3); %OxyRn
f0=param(4); %OxyRa
g0=param(5); %Pkat
h0=param(6); %mRNAtev
j0=param(7); %TEVp1
k0=param(8); %GFPt - alive
l0=param(9); %GFPd
z0=param(10); %RFPt
x0=param(11); %RFPd - alive
N0 = [a0;s0;d0;f0;g0;h0;j0;k0;l0;z0;x0];

% Rate Constants (Mol/sec)
k1 = 0.00009627; %cannot change
k2 = 0.00042; %cannot change
k3 = 0.00001;
k4 = 0.00002;
k5 = 0.000001;
k6 = 0.0008;
alpha = 0.7; 
k7 = 0.0001; 

m1 = 0;
m2 = 0.01;
m3 = 0.01;
m4 = 0.08;

d1 = 0.00000802; %cannot change
d2 = 0.001;
d3 = 0.001;
d4 = 0.009;
d5 = 0.0021; %cannot change
d6 = 0.0000688; %cannot change
d7 = 0.00000963; %cannot change
d8 = 0.0000642; %cannot change
d9 = 0.000138; %cannot change

%Running the Built-in ODE Solver
[tt,NN] = ode45(@f,[0,tf],N0);

AA = NN(:,1); %Diamine
SS = NN(:,2);
DD = NN(:,3);
FF = NN(:,4);
GG = NN(:,5);
HH = NN(:,6);
JJ = NN(:,7);
KK = NN(:,8);% GFP alive
LL = NN(:,9); 
ZZ = NN(:,10);
XX = NN(:,11);% RFP alive

%figure()
%plot(tt,AA,tt,SS,tt,DD,tt,FF,tt,GG,tt,HH, tt,JJ, tt, KK, tt, LL, tt, ZZ, tt, XX);
%legend('Diamine','H2O2', 'OxyRn', 'OxyRa', 'Pkat', 'mRNAtev', 'TEVp1', 'GFPt', 'GFPd', 'RFPt', 'RFPd');
%title('Everything');
%xlabel('BA');
%ylabel('Concentration');
%ylim([0 50]);

%figure()
%plot(tt,AA,tt,SS,tt,DD,tt,FF);
%legend('Diamine','H2O2', 'OxyRn', 'OxyRa');
%title('Part 1');
%xlabel('BA');
%ylabel('Concentration');
%ylim([0 200]);

%figure()
%plot(tt,GG,tt,HH, tt,JJ);
%legend('Pkat', 'mRNAtev', 'TEVp1');
%title('Part 2');
%xlabel('BA');
%ylabel('Concentration');
%ylim([0 2]);

%figure()
%plot(tt, KK, 'green', tt, LL, tt, ZZ, tt, XX, 'red');
%legend('GFPt', 'GFPd', 'RFPt', 'RFPd');
%title('Part 3');
%xline(500000);
%xlabel('BA');
%ylabel('Concentration');
%ylim([0 120]);

%TARGET GRAPH
figure()
plot(tt, KK, 'green', tt, XX, 'red');
%legend('GFPt', 'RFPd');
title('Target Graph');
xline(1000, '-','Point of Inspection');
xlabel('Time');
ylabel('Fluorescence Concentration');
ylim([0 30]);
xlim([0 1500]);

    function dYdt = f(t,Y)

        Dia=Y(1);
        H2O2=Y(2);
        OxyRn=Y(3);
        OxyRa=Y(4);
        Pkat=Y(5);
        mRNAtevp=Y(6);
        TEVp1=Y(7);
        GFPt=Y(8);
        GFPd=Y(9);
        RFPt=Y(10);
        RFPd=Y(11);
        
        %ODEs that govern the concentrations
        dDiadt = m1 - k1.*Dia;
        dH2O2dt = k1.*Dia - k2.*OxyRn.*H2O2 - d1.*H2O2;
        dOxyRndt = m2 - k2.*OxyRn.*H2O2 - d2.*OxyRn;
        dOxyRadt = k2.*OxyRn.*H2O2 - d3.*OxyRa;
	    dPkatdt = - k3.*OxyRa.*Pkat;
        dmRNAtevpdt = alpha.*(OxyRa./k7)./(1+(OxyRa./k7))- d4.*mRNAtevp;


	    dTEVp1dt = k4.*mRNAtevp - d7.*TEVp1;
	    dGFPtdt = m3 - k5.*GFPt.*TEVp1 - d5.*GFPt; %!!
	    dRFPtdt = m4 - k6.*RFPt.*TEVp1 - d9.*RFPt;
        dGFPddt = k5.*GFPt.*TEVp1 - d8.*GFPd;
        dRFPddt = k6.*RFPt.*TEVp1 - d6.*RFPd; %!!
      

        dYdt=[dDiadt;dH2O2dt;dOxyRndt;dOxyRadt;dPkatdt;dmRNAtevpdt;dTEVp1dt;dGFPtdt;dGFPddt;dRFPtdt;dRFPddt];
        

    end


end