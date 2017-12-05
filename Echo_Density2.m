function [theta] = Echo_Density2(h,t,fs)

% load ('~/Desktop/reverberation/MMR ir data/mmr/mmr1.mat')
% % [h,fs] = audioread('mmr1.wav');
% h = y';

t1 = 1000; % start point 
t2 = length(h); % total time
% t2 = 0.1*fs; % total time
L = round(t*fs/2); % half window length
w = hann(L*2+1); % window function
w = w/sum(w); % normalize the window
theta = zeros(1,t2); % profile
choice = 1;

for n = 1:1:t2-2*L-1
    if choice == 1
        segma = sqrt(sum(w.*h(n:n+2*L).^2));
        theta(n) = sum(w.*(abs(h(n:n+2*L))>segma));
    else
        segma = sqrt(sum(h(n:n+2*L).^2)/(2*L+1));
        theta(n) = sum((abs(h(n:n+2*L))>segma))/(2*L+1);
    end
end
theta = theta/erfc(1/sqrt(2));
theta = [zeros(1,L) theta];

% figure
% t = t1:t2;
% plot(t,h(t)*50);
% hold on
% plot(t,theta(t));
% hold off
% set(gca,'XScale','log')
% grid;

