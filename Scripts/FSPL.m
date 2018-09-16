function AttenuationFS = FSPL(distance,freq)
%Calculating the Free Space Path Loss for a HF wave.
AttenuationFS = 20*log10(freq) + 20*log10(distance) + 32.45;
end

