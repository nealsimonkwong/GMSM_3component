function [Y, sig_lnY] = gmpmVoverH_BC2016(M,Rrup,Rjb,Rx,lambda,dip,Fhw,Vs30,region,Sj,T,W,Ztor,Z2p5,Zhyp)
% By: Neal Simon Kwong, nealsimonkwong@berkeley.edu
% Last modified: 4/27/18
% 
% Bozorgnia and Campbell 2016 GMPM for V/H ratios. Paper citation:
% 
% Bozorgnia, Y. and Campbell, K. W. (2016). Ground Motion Model 
% for the Vertical-to-Horizontal (V/H) Ratios of PGA, PGV, and Response
% Spectra. Earthquake Spectra, 32(2), 951-978.
% 
% Required inputs:
% M       = Moment magnitude
% Rrup    = Closest distance to rupture plane (km)
% Rjb     = Closest distance to surface projection of rupture plane (km)
% Rx      = Horizontal distance from surface projection of top edge of 
%         rupture plane to site, measured perpendicular to average strike, 
%         and is negative for footwall but positive for hanging wall (km)
% lambda  = rake angle (degree); average angle of slip measured in
%         the plane of rupture; used to compute FRV and FNM
% dip     = Average dip angle of rupture plane (deg)
% Fhw     = hanging wall effect (=1 when included and =0 when excluded)
% Vs30    = Time-averaged shear wave velocity in top 30m of site (m/s)
% region  = 0 for Global
%         = 1 for California
%         = 2 for Japan and Italy
%         = 3 for eastern China
% Sj      = Indicator variable representing Japan's site effects
% T       = Vibration period of interest (sec); =0 for PGA and =-1 for PGV
% 
% Optional inputs (default values can be calculated):
% W       = Down-dip width of rupture plane (km)
% Ztor    = Depth to top of rupture plane (km)
% Z2p5    = Depth to 2.5km/s shear wave velocity, or sediment depth (km)
% Zhyp    = Hypocentral depth of earthquake measured from sea level (km)
% 
% Outputs:
% Y       = Median V/H ratio for: 5%-damped spectral acceleration (g) or PGA
%         (g) when T=0 or PGV (cm/s) when T=-1
% sig_lnY = Total aleatory standard deviation
% 
% Other variables:
% FRV     = Indicator variable representing reverse and reverse-oblique
%         faulting (=1 when rake is within 30 to 150 deg)
% FNM     = Indicator variable representing normal and normal-oblique
%         faulting (=1 when rake is within -30 to -150 deg)
% Zbot    = Depth to the bottom of the seismogenic crust (km); needed only
%          when W is unknown

%% If desired, can check inputs for validity here (e.g., Rrup>0)

%% Calculate remaining input variables for all GMPMs
% Style of faulting
FRV = (lambda > 30 & lambda < 150);
FNM = (lambda > -150 & lambda < -30);
% Depth to bottom of seismogenic crust (km); needed for GMPM for H comp. when W unknown
Zbot = 15; 
% Calculate default values for input vars (see Excel file)
if nargin<12
    Zbor = 15; % Depth to bottom of rupture plane (km)
    Ztor = FRV*max([2.704-1.226*max([M-5.849 0]) 0])^2 + (1-FRV)*max([2.673-1.136*max([M-4.97 0]) 0])^2 ;
    if Ztor>20; display('Warning: Ztor exceeds 20km.'); end
    W = min([sqrt(10^((M-4.07)/0.98)) (Zbot-Ztor)/sind(dip)]);
    Z2p5 = 0.6068;
    temp1 = (M<6.75)*(-4.317+0.984*M) + (M>=6.75)*2.325 + (dip<40)*(0.0445*(dip-40));
    temp2 = log(0.9*(Zbor-Ztor));
    Zhyp = max([Ztor+exp(min([temp1 temp2])) 5]);
    if Zhyp>20; display('Warning: Zhyp exceeds 20km.'); end
end

%% Call external GMPMs for input periods
% Use GMPM for H comp. (from Prof. Baker)
[Y_H, ~, ~, tau_lnY_H, phi_lnY_H] = CB_2014_nga(M, T, Rrup, Rjb, Rx, W, Ztor, Zbot, dip, lambda, Fhw, Vs30, Z2p5, Zhyp, region);
% Use GMPM for V comp.
[Y_V, ~, tau_lnY_V, phi_lnY_V] = gmpmV_BC2016(M,Rrup,Rjb,Rx,FRV,FNM,dip,Vs30,region,Sj,T,W,Ztor,Z2p5,Zhyp);
% Use correlation model for H and V components at same vibration period
[rhoB_V_H, rhoW_V_H] = corrModel_HVsameT_BC2016(T,M);

%% Ensure outputs from all GMPMs are same size
Y_H = reshape(Y_H,size(Y_V));
tau_lnY_H = reshape(tau_lnY_H,size(Y_V));
phi_lnY_H = reshape(phi_lnY_H,size(Y_V));

%% Median prediction
Y = Y_V./Y_H;

%% Aleatory variability
tau_lnY = sqrt( tau_lnY_V.^2 + tau_lnY_H.^2 - 2*rhoB_V_H .* tau_lnY_V .* tau_lnY_H );
phi_lnY = sqrt( phi_lnY_V.^2 + phi_lnY_H.^2 - 2*rhoW_V_H .* phi_lnY_V .* phi_lnY_H );
sig_lnY = sqrt( tau_lnY.^2 + phi_lnY.^2 );

end