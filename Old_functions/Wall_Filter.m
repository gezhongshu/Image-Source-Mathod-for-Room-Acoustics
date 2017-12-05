function N = Wall_Filter(octave_band, walltype, fs, DispFilter) 
walltype = sqrt(1-walltype);
octave_band = octave_band/(fs/2);
curve = fit( octave_band', walltype', 'exp2');
cutoff = 1;
for n = 0:0.02:1
    if n > 8000/(fs/2)
        cutoff = n-0.02;
        scale = curve(cutoff);
        break
    end
end

N = 20;
B1 = 0:0.02:cutoff;
A1 = curve(B1)';
B2 = cutoff+0.02:0.02:1;
A2 = ones(size(B2))*scale;

F = [B1 B2];
A = [A1 A2];
d = fdesign.arbmag('N,F,A',N,F,A);
filter = design(d,'freqsamp','SystemObject',true);
% extract the numerator of the filter
s = coeffs(filter);
N = s.Numerator;
if DispFilter==1
fvtool(filter,'MagnitudeDisplay','Zero-phase','Color','White');
end