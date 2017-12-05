t = 0:.01:pi;
t2 = 0:-0.01:-pi;
polar(t,1-(1-0.1)*t/pi,'--b');
hold on;
polar(t2,1+(1-0.1)*t2/pi,'--b');
hold off