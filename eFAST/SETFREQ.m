% Algorithm for selection of a frequency
% set for the complementary group. Done
% recursively as described in:
% Appendix of "Sensitivity Analysis"
% [Saltelli et al., 2000]
%OMci = SETFREQ(k-1,OMi/2/MI) 
function OMci = SETFREQ(Kci,OMciMAX,f)
% Kci = k = 23
% OMciMAX = 0.6000
% f = i = 1-23

if Kci==1
    OMci = 1;
elseif OMciMAX==1
    OMci = ones(1,Kci);
else
    if(OMciMAX < Kci)
         INFD = OMciMAX;
         % INFD = 0.6000
    else
        INFD = Kci;
    end
    
    ISTEP = round((OMciMAX-1)/(INFD-1));
    % ISTEP = 1
    
    if(OMciMAX == 1)
        ISTEP = 0;
    end
    % Must be deleted.
    OTMP = 1:ISTEP:INFD*ISTEP;
    % OTMP = 1:1:0.6;

    fl_INFD = floor(INFD);
    for i=1:Kci
        j = mod(i-1,fl_INFD)+1;
        OMci(i) = OTMP(j);
    end
end
OMci(f)=[];