function [rhoB_V_H, rhoW_V_H] = corrModel_HVsameT_BC2016(T,M)
% By: Neal Simon Kwong, nealsimonkwong@berkeley.edu
% Last modified: 4/27/18
% 
% Bozorgnia and Campbell 2016 GMPM for correlations between vertical and
% horizontal epsilons at the same vibration period. Paper citation:
% 
% Bozorgnia, Y. and Campbell, K. W. (2016). Ground Motion Model 
% for the Vertical-to-Horizontal (V/H) Ratios of PGA, PGV, and Response
% Spectra. Earthquake Spectra, 32(2), 951-978.
% 
% Required inputs:
% T       = Vibration period of interest (sec); =0 for PGA and =-1 for PGV
% M       = Moment magnitude
%
% Outputs:
% rhoB_V_H  = Between-event correlation 
% rhoW_V_H  = Within-event correlation 

%% Specify periods used in GMPM development
T_BC2016 = [0.010 0.020 0.030 0.050 0.075 0.10 0.15 0.20 0.25 0.30 0.40 0.50 0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0 0 -1]; % in units of sec except =0 for PGA and =-1 for PGV
T_BC2016sub = T_BC2016(1:(end-2)); % Remove periods for PGA and PGV

%% Execute GMPM for input periods
nTinput = length(T);
rhoB_V_H = zeros(size(T));
rhoW_V_H = zeros(size(T));
for ii=1:nTinput
    Tcurr = T(ii);
    if ismember(Tcurr,T_BC2016) % Pre-defined period requested
        % Report prediction for pre-defined period
        ip = find(T_BC2016 == Tcurr);
        [rhoB_V_H(ii), rhoW_V_H(ii)] = BC2016corr_sub(ip,M);
    elseif Tcurr>=min(T_BC2016sub)&& Tcurr<=max(T_BC2016sub) % Assume Tcurr is between 0.01 and 10 sec
        % Determine neighboring periods        
        ip_lo = find(T_BC2016sub<=Tcurr,1,'last');
        ip_hi = find(T_BC2016sub>=Tcurr,1,'first');
        T_lo = T_BC2016sub(ip_lo);
        T_hi = T_BC2016sub(ip_hi);
        
        % Determine output for neighboring periods
        [rhoB_V_H_lo, rhoW_V_H_lo] = BC2016corr_sub(ip_lo,M);
        [rhoB_V_H_hi, rhoW_V_H_hi] = BC2016corr_sub(ip_hi,M);
        
        % Linearly interpolate
        rhoB_V_H(ii) = interp1(log([T_lo T_hi]),[rhoB_V_H_lo rhoB_V_H_hi],log(Tcurr)); % Semilog scale (can also try log-log scale by taking log of rho's)
        rhoW_V_H(ii) = interp1(log([T_lo T_hi]),[rhoW_V_H_lo rhoW_V_H_hi],log(Tcurr)); % Semilog scale (can also try log-log scale by taking log of rho's)        
    else
        display(['Current input period is ' num2str(Tcurr,4)]);
        error('Invalid input period');
    end
end
end


function [rhoB_V_H, rhoW_V_H] = BC2016corr_sub(ip,M)
%% Regression output
% Col headers: rhoW1 rhoW2 rhoB1 rhoB2
Table1 = [...
    0.783	0.718	0.916	0.893
    0.785	0.718	0.917	0.891
    0.784	0.703	0.919	0.884
    0.782	0.678	0.931	0.877
    0.781	0.681	0.933	0.909
    0.778	0.657	0.944	0.9
    0.775	0.653	0.946	0.875
    0.774	0.63	0.946	0.826
    0.77	0.642	0.945	0.784
    0.757	0.658	0.949	0.819
    0.75	0.661	0.953	0.719
    0.742	0.666	0.955	0.716
    0.73	0.688	0.965	0.718
    0.721	0.684	0.967	0.765
    0.707	0.663	0.967	0.796
    0.691	0.664	0.966	0.799
    0.652	0.629	0.956	0.822
    0.669	0.613	0.945	0.839
    0.665	0.586	0.921	0.86
    0.586	0.573	0.898	0.685
    0.639	0.547	0.854	0.72
    0.782	0.72	0.915	0.893
    0.754	0.68	0.882	0.699];

%% Between-event
idCol_B = 3;
if M<=4.5
    rhoB_V_H = Table1(ip,idCol_B);
elseif M>=5.5
    rhoB_V_H = Table1(ip,idCol_B+1);    
else
    val1 = Table1(ip,idCol_B);
    val2 = Table1(ip,idCol_B+1);  
    rhoB_V_H = val2 + (val1-val2)*(5.5-M);
end

%% Within-event
idCol_W = 1;
if M<=4.5
    rhoW_V_H = Table1(ip,idCol_W);
elseif M>=5.5
    rhoW_V_H = Table1(ip,idCol_W+1);    
else
    val1 = Table1(ip,idCol_W);
    val2 = Table1(ip,idCol_W+1);  
    rhoW_V_H = val2 + (val1-val2)*(5.5-M);
end

end