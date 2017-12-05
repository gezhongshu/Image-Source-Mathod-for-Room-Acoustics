load ('~/Desktop/reverberation/room-measures/source.mat')
load('~/Desktop/reverberation/room-measures/mmr/mmr1response.mat')
addpath /Users/silver/Documents/MATLAB/Wen_ISM
y = fftdeconv(response, x, 0.99);
y = y(1000:fs)';

echo = Echo_Density(y,fs);
echo2 = Echo_Density2(y,0.01,fs);
echo3 = Echo_Density2(y,0.03,fs);
echo(1:537) = 0;
echo2(1:735) = 0;
% echo3(1:537) = 0;

figure;
T = 1:size(y,1);
semilogx(T,y*35)
hold on
semilogx(T,echo2(T))
semilogx(T,echo(T))
semilogx(T,echo3(T))
xlabel('Time (ms)')
ylabel('Amplitude/Echo Density')
axis([500 45000 -0.4 1.4 ]);
legend('Measured Impulse Response','Echo Density Profile Which Window Size is 0.01s','Echo Density Profile Which Window Size is 0.02s', 'Echo Density Profile Which Window Size is 0.03s')
grid;
hold off