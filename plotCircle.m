%%Plotting a circle to represent Venus. Helps give perspective to
%%trajectories

x(1) = 0;
y(1) = 0;

theta = 0:1:360;
y = sind(theta)*R;
x = cosd(theta)*R;
plot(x,y,'r');
legend('Venus');