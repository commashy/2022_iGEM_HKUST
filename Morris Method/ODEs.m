 function answer = ODEs(param)
% reaction_12(30000,[0,0,0,0,0,0])
t = 30000;
% initial concentrations
a0=0.01; %Diamine
s0=0; %H2O2
d0=0; %OxyRn
f0=0; %OxyRa
g0=1; %Pkat
h0=1; %PoxyS
j0=0; %GFP
k0=0; %RFP
z0=0; %mRNAGFP
x0=0; %mRNARFP

tf = t.*2;

k1 = param(1); 
k2 = param(2); 
k3 = param(3);
k4 = param(4);
k5 = param(5);
k6 = param(6);
k7 = param(7);
alpha = param(8); 

m1 = param(9);
m2 = param(10);
m3 = param(11);
m4 = param(12);

d1 = param(13);
d2 = param(14);
d3 = param(15);
d4 = param(16); 
d5 = param(17); 
d6 = param(18); 
d7 = param(19); 
d8 = param(20); 

%rate constants, 1/sec 1/mM

% initial concentration
N0 = [a0;s0;d0;f0;g0;h0;j0;k0;z0;x0];

%opt=odeset('Events',@ReachSS);
%option2 = odeset('NonNegative',1);
[tt,NN] = ode15s(@f,[0,tf],N0);

AA = NN(:,1);
SS = NN(:,2);
DD = NN(:,3);
FF = NN(:,4);
GG = NN(:,5);
HH = NN(:,6);
JJ = NN(:,7); %GFP
KK = NN(:,8); %RFP
ZZ = NN(:,9); 
XX = NN(:,10);
answer = KK(end);

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
        dDiadt = m1 - k1.*Dia;
        dH2O2dt = k1.*Dia - k2.*OxyRn.*H2O2 - d1.*H2O2;
        dOxyRndt = m2 - k2.*OxyRn.*H2O2 - d2.*OxyRn;
        dOxyRadt = k2.*OxyRn.*H2O2 - d3.*OxyRa;
	    dPkatdt = - k3.*OxyRa.*Pkat;
        dmRNAtevpdt = alpha.*(OxyRa./k7)./(1+(OxyRa./k7)) - k4.*mRNAtevp - d4.*mRNAtevp;

	    dTEVp1dt = k4.*mRNAtevp - d7.*TEVp1;
	    dGFPtdt = m3 - k5.*GFPt.*TEVp1 - d5.*GFPt; %!!
	    dRFPtdt = m4 - k6.*RFPt.*TEVp1 - d9.*RFPt;
        dGFPddt = k5.*GFPt.*TEVp1 - d8.*GFPd;
        dRFPddt = k6.*RFPt.*TEVp1 - d6.*RFPd; 

        dYdt=[dDiadt;dH2O2dt;dOxyRndt;dOxyRadt;dPkatdt;dPoxySdt;dGFPdt;dRFPdt;dmRNAGFPdt;dmRNARFPdt];


    end

    function [position, isterminal, direction] = ReachSS(t,N)
        position = max(abs(f(t,N))<0.001) - 0.5;
        isterminal = 1;
        direction = 0;
    end

end
