function polarization = polarization(hOrV,sigma,permittivity,conductivity,freq)
wavelength = 3*10^8 / (freq * 10^6);
permPrime = permittivity - 60* wavelength * conductivity * 1i;
srRoot = sqrt(permPrime - (cos(sigma))^2);
%Returns the polarization of the HF wave.
if hOrV == 'h'
    polarization = (sin(sigma) - srRoot)/(sin(sigma) + srRoot);
elseif hOrV == 'v'
    polarization = (permPrime* sin(sigma) - srRoot)/(permPrime*sin(sigma) + srRoot);
end

