addpath /Users/silver/Documents/MATLAB/Wen_ISM
load('/Users/silver/Documents/MATLAB/Wen_Ray_Tracing/RT_MMR_1KHz_Pos1_270*270');
dt = 0.01; % time interval in energy response (in second)
fs = 48000;
T =  0:dt:0.1-dt;
E =  Energy_resp(IR, dt, fs);

figure;
stem(T, E, 'LineWidth',3,'Color',[0 0 0]);
xlabel('Time(s)')
ylabel('Sum of the energy')
axis([0 0.1 0 25000]);
grid;