
clc, clear;

[h,fs] = audioread('mmr1.wav');
t1 = 1000; % start point 
t2 = 1*fs; % total time
L = round(0.02*fs/2); % half window length
w = hann(L*2+1); % window function
w = w/sum(w); % normalize the window
theta = zeros(1,t2); % profile
choice = 1;

for n = t1:1:t2
    if choice == 1
        segma = sqrt(sum(w.*h(n-L:n+L).^2));
        theta(n) = sum(w.*(abs(h(n-L:n+L))>segma));
    else
        segma = sqrt(sum(h(n-L:n+L).^2)/(2*L+1));
        theta(n) = sum((abs(h(n-L:n+L))>segma))/(2*L+1);
    end
end
theta = theta/erfc(1/sqrt(2));

figure
t = t1:t2;
plot(t,h(t)*50);
hold on
plot(t,theta(t));
hold off
set(gca,'XScale','log')
grid;

