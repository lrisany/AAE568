function projecttpbvp3_ext()
clc 
clear all

t0 = 0;
tf = 1;

solinit = bvpinit(linspace(t0, tf), [2;2;2;2;2;2;.25;.25;.25;.25;.25;.25;2;2;2]);
options = bvpset('Stats', 'on', 'RelTol', 1e-1);
solution = bvp4c(@BVP_ode_project3, @boundary_conditions_project3, solinit, options);

I1 = 396.2;
I2 = 1867;
I3 = 1987.8;

y = solution.y;
t = solution.x;


B1 = y(1,:);
B2 = y(2,:);
B3 = y(3,:);
B4 = y(4,:);

Gamma1 = y(5,:);
Gamma2 = y(6,:);
Gamma3 = y(7,:);
Gamma4 = y(8,:);


w1 = y(9,:);
w2 = y(10,:);
w3 = y(11,:);

lambda1 = y(12,:);
lambda2 = y(13,:);
lambda3 = y(14,:);

L1 = -lambda1 / I1;
L2 = -lambda2 / I2;
L3 = -lambda3 / I3;

hold on
plot(t, w1, 'b');
plot(t, w2, 'r');
plot(t, w3, 'g');
title('Angular Velocity Time Histories');
xlabel('time (s)');
ylabel('omega (rad/s)');
legend('w1', 'w2', 'w3');
grid on

figure
hold on

plot(t, lambda1, 'b--');
plot(t, lambda2, 'r--');
plot(t, lambda3, 'g--');
title('Costate Time Histories');
xlabel('time (s)');
ylabel('Costates');
legend('Lambda 1', 'Lambda 2', 'Lambda 3');

figure;
hold on
plot(t, L1, 'b');
plot(t, L2, 'r');
plot(t, L3, 'g');
title('Time Histories');
xlabel('time (s)');
ylabel('Control Torques (Nm)');
legend('L1', 'L2', 'L3');
grid on

figure;
hold on
plot(t, B1, 'b');
plot(t, B2, 'r');
plot(t, B3, 'g');
plot(t, B4, 'k');
title('Euler Parameters');
xlabel('time (s)');
ylabel('Euler Paramaters');
legend('E1', 'E2', 'E3', 'E4');
grid on

hold off