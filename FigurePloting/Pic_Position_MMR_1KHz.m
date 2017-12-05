fs = 48000;
load('/Users/silver/Documents/MATLAB/Wen_ISM/ISM_MMR_1KHz_Pos1.mat');
ISM_1 = IR';

load('IR_MMR_1KHz_Pos1.mat');
Measure_1 = y;

figure;
hold on
plot((0:size(Measure_1,2)-1)/fs,Measure_1(1:size(Measure_1,2)))
plot((0:size(ISM_1,2)-1)/fs,ISM_1(1:size(ISM_1,2)),'LineWidth',1)
xlabel('Time (s)')
ylabel('Power')
title('Position 1')
axis([0 0.12 -0.00015 0.00015]);
grid;
legend('Measured Impulse Response','The Image Sourse Model')
hold off

echo_Measure_1 = Echo_Density(Measure_1',fs);
echo_Measure_1(1:537) = 0;
echo_ISM_1 = Echo_Density(ISM_1',fs);
figure;
hold on
plot((0:size(echo_Measure_1,2)-1)/fs, echo_Measure_1)
plot((0:size(echo_ISM_1,2)-1)/fs, echo_ISM_1)
xlabel('Time (s)')
ylabel('Echo Density')
title('Position 1')
axis([0 0.11 0 1.2]);
legend('Measured Impulse Response','The Image Sourse Model')
hold off

% Position 2--------------------------------------------------------------%
load('/Users/silver/Documents/MATLAB/Wen_ISM/ISM_MMR_1KHz_Pos2.mat');
ISM_2 = IR';

load('IR_MMR_1KHz_Pos2.mat');
Measure_2 = y;

figure;
hold on
plot((0:size(Measure_2,2)-1)/fs,Measure_2(1:size(Measure_2,2)))
plot((0:size(ISM_2,2)-1)/fs,ISM_2(1:size(ISM_2,2)),'LineWidth',1)
xlabel('Time (s)')
ylabel('Power')
title('Position 2')
axis([0 0.12 -0.00015 0.00015]);
grid;
legend('Measured Impulse Response','The Image Sourse Model')
hold off

echo_Measure_2 = Echo_Density(Measure_2',fs);
echo_Measure_2(1:885) = 0;
echo_ISM_2 = Echo_Density(ISM_2',fs);
figure;
hold on
plot((0:size(echo_Measure_2,2)-1)/fs, echo_Measure_2)
plot((0:size(echo_ISM_2,2)-1)/fs, echo_ISM_2)
xlabel('Time (s)')
ylabel('Echo Density')
title('Position 2')
axis([0 0.11 0 1.2]);
legend('Measured Impulse Response','The Image Sourse Model')
hold off

% Position 3--------------------------------------------------------------%
load('/Users/silver/Documents/MATLAB/Wen_ISM/ISM_MMR_1KHz_Pos3.mat');
ISM_3 = IR';

load('IR_MMR_1KHz_Pos3.mat');
Measure_3 = y;

figure;
hold on
plot((0:size(Measure_3,2)-1)/fs,Measure_3(1:size(Measure_3,2)))
plot((0:size(ISM_3,2)-1)/fs,ISM_3(1:size(ISM_3,2)),'LineWidth',1)
xlabel('Time (s)')
ylabel('Power')
title('Position 3')
axis([0 0.12 -0.00015 0.00015]);
grid;
legend('Measured Impulse Response','The Image Sourse Model')
hold off

echo_Measure_3 = Echo_Density(Measure_3',fs);
echo_Measure_3(1:748) = 0;
echo_ISM_3 = Echo_Density(ISM_3',fs);
figure;
hold on
plot((0:size(echo_Measure_3,2)-1)/fs, echo_Measure_3)
plot((0:size(echo_ISM_3,2)-1)/fs, echo_ISM_3)
xlabel('Time (s)')
ylabel('Echo Density')
title('Position 3')
axis([0 0.11 0 1.2]);
legend('Measured Impulse Response','The Image Sourse Model')
hold off