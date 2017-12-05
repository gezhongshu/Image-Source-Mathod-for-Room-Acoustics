function [B,A] = BPFilter(Fc, Fs, N ,DispFilter) 

% f1 = Fc/(2^(1/6)); 
% f2 = Fc*(2^(1/6)); 
% Qr = Fc/(f2-f1); 
% Qd = (pi/2/N)/(sin(pi/2/N))*Qr;
% alpha = (1 + sqrt(1+4*Qd^2))/2/Qd; 
% W1 = Fc/(Fs/2)/alpha; 
% W2 = Fc/(Fs/2)*alpha;
% [B,A] = butter(N,[W1,W2],'bandpass'); 

f1 = Fc/2; 
f2 = Fc*(2^(1/6)); 
Wp = [f1 f2]/(Fs/2);
% [B,A] = butter(N,Wp,'bandpass');
% [B,A] = ellip(N,5,40,Wp);
[B,A] = cheby1(N,0.2,Wp);


if DispFilter==1
freqz(B,A,Fs);
end