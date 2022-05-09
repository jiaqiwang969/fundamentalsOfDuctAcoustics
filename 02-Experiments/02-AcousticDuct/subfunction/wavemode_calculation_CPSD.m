%i_method = 1: Ê±¼äÆ×
%i_method = 2: ½×±ÈÆ×
function [GAMMA,freq]=wavemode_calculation_CPSD(data,fs,nk)
    signal = data(:,1:nk);
    
        
    L_signal = length(signal);
    L_seg = round(L_signal/10);
    Wind = hamming(L_seg);
    Noverlap = round(L_seg/2);
    Nfft = 2^(ceil(log2(L_seg))+1);  
    for k=1:nk
        for l = 1:nk
            [C{k,l},freq] = cpsd(signal(:,k),signal(:,l),Wind,Noverlap,Nfft,fs);          
        end
    end
    GAMMA = zeros(Nfft/2+1,nk+1);
    mode=-nk/2:nk/2;
    for m = 1:nk+1
        temp_f = zeros(Nfft/2+1,1);
        for k = 1:nk
            for l = 1:nk
                temp_f = temp_f + 0.5*C{k,l}*exp(i*mode(m)*2*pi*k/nk)*exp(-i*mode(m)*2*pi*l/nk);
            end
        end
        GAMMA(:,m) = temp_f/(nk*nk);
    end
    
    GAMMA = 10*log10(abs(GAMMA)/4e-10);
%   GAMMA = abs(GAMMA);
    h=figure
    imagesc([-32/2:32/2],freq,GAMMA); 
 end
   
