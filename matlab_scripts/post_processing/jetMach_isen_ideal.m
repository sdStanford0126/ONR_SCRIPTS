function Mj = jetMach_isen_ideal(NPR,gamma)
% calculate ideal mach jet number assuming isentropic ideal expansion based
% on NPR
% input: NPR: nozzle pressure ratio (total pressure to ambient)
%        gamma: specific heat ratio
% output: jet Mach number
% ref: https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
%Mj =  sqrt(2/(gamma-1)*((NPR)^((gamma-1)/gamma)-1));
% actuaaly solving for acoustic mach number
Mj = sqrt(2/(gamma-1)*(1-(1/NPR)^((gamma-1)/gamma)));