fs = 48000;
load('/Users/silver/Documents/MATLAB/Wen_ISM/ISM_Tanna_1KHz_Pos1.mat');
ISM = IR';
load('/Users/silver/Documents/MATLAB/Wen_Ray_Tracing/RT_Tanna_1KHz_Pos1_270*270');
RT = IR';
% load('/Users/silver/Documents/MATLAB/Wen_Ray_Tracing/RT_Tanna_1KHz_Pos1_270*270_3rd');
% RT_3rd = IR';
% load('/Users/silver/Documents/MATLAB/Wen_Ray_Tracing/RT_Tanna_1KHz_Pos1_270*270_4th');
% RT_4th = IR';

load('IR_Tanna_1KHz_Pos1.mat');
Measure = y;

[ma,R]=max(Measure); % Find the direct soudn power, which is the lasrgest one
RT=RT*ma/max(RT); % normalize the model energy and the measured energy
% RT_3rd=RT_3rd*ma/max(RT_3rd); % normalize the model energy and the measured energy
% RT_4th=RT_4th*ma/max(RT_4th); % normalize the model energy and the measured energy

figure;
hold on
plot((0:size(Measure,2)-1)/fs,Measure(1:size(Measure,2)))
plot((0:size(ISM,2)-1)/fs,ISM(1:size(ISM,2)),'LineWidth',1)
plot((0:size(RT,2)-1)/fs,RT(1:size(RT,2)),'LineWidth',1)
xlabel('Time (s)')
ylabel('Power')
axis([0 0.12 -0.0002 0.0002]);
grid;
legend('Measured Impulse Response','The Image Sourse Model','The Ray-Tracing Method')
hold off

echo_Measure = Echo_Density(Measure',fs);
echo_Measure(1:12) = 0;
echo_ISM = Echo_Density(ISM',fs);
echo_RT = Echo_Density(RT',fs);
% echo_RT_3rd = Echo_Density(RT_3rd',fs);
% echo_RT_4th = Echo_Density(RT_4th',fs);

figure;
hold on
plot((0:size(echo_Measure,2)-1)/fs, echo_Measure)
plot((0:size(echo_ISM,2)-1)/fs, echo_ISM)
plot((0:size(echo_RT,2)-1)/fs, echo_RT)
% plot((0:size(echo_RT_3rd,2)-1)/fs, echo_RT_3rd)
% plot((0:size(echo_RT_4th,2)-1)/fs, echo_RT_4th)
axis([0 0.11 0 1.2]);
% set(gca,'xtick',[0 1 2 24]);
% set(gca,'ytick',[0 100 200 300 400 500 600 700 800 900]);
legend('Measured Impulse Response','The Image Source Model','The Ray-Tracing Method')
hold off