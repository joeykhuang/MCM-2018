SNRmodes = ['f', 's', 'n';'r', 's', 'n';'r', 's', 'y';'f', 'g', 'n';'r', 'g', 'y';'r', 'g', 'n';];
seaParams = [80, 4];
landParams = [13, 0.004];

for i=1:6
filename = [SNRmodes(i,1),SNRmodes(i,2),SNRmodes(i,3),'Prime.csv'];
SNR = [0 0 0 0 0 0];
    for j=1:6
        if SNRmodes(i,2) == 's'
        SNR(j) = SNR(j) + AttenuationTotal(j,'E',1000,seaParams(1),seaParams(2),5,100,SNRmodes(i,1),SNRmodes(i,2),SNRmodes(1,3));
        elseif SNRmodes(i,2) == 'g'
        SNR(j) = SNR(j) + AttenuationTotal(j,'E',1000,landParams(1),landParams(2),5,100,SNRmodes(i,1),SNRmodes(i,2),SNRmodes(1,3));
        end
    end
csvwrite(filename,SNR)
end