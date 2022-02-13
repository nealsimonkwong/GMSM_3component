%% Calculate pseudo-acceleration response spectrum:
% Damping is specified as input argument
% Method in Prof Chopra's book is encoded as "LinearInterpolation"
function RS = calcResponseSpectrum(ugdd,dt,T,zeta)
RS = zeros(length(T),1); % RS = response spectrum
for ii = 1:length(T)
    clear disp_history
    wn = (2*pi)./T(ii);
    disp_history = LinearInterpolation(1,wn,zeta,0,0,ugdd,dt,1);
    RS(ii,1) = max(abs(disp_history))*wn^2;
end
end

