fs = 48000;
sourceName = 'source.mat';
responseName = 'pollack/pollack1a.mat';
N = 2^20;

load( sourceName );
load( responseName );

y = fftdeconv(response, x, 0.99);

maxTime = 3.0;
maxN = round( maxTime * fs );
y = y( 1:maxN ); % cut result to maxTime length

t = (0:length(y)-1) / fs;
plot(t, y)
xlabel('Time (seconds)');
title('Impulse Response');
