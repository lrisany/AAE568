clear; clc; close all;

% Initial Conditions
w1_i = 1;
w2_i = -0.1;
w3_i = 0.1;
e1_i = 0;     % Euler Parameters
e2_i = 0;
e3_i = 0;
e4_i = 1;

% Initial Condition inputs
init = [w1_i w2_i w3_i e1_i e2_i e3_i e4_i];
tSpan = 0:0.01:30;

% Numerical Integration
value = 1e-6;
options = odeset('reltol', value, 'abstol', value);
[tAll, f] = ode45(@(t,x) eom(x), tSpan, init, options);

% Parse outputs
omega1All = f(:,1);
omega2All = f(:,2);
omega3All = f(:,3);

e1 = f(:,4);
e2 = f(:,5);
e3 = f(:,6);
e4 = f(:,7);

C11 = 1 - 2*e2.^2 - 2*e3.^2;
C12 = 2 * (e1.*e3 - e3.*e4);
C13 = 2 * (e3.*e1 + e2.*e4);
C21 = 2 * (e1.*e2 + e3.*e4);
C22 = 1 - 2*e3.^2 - 2*e1.^2;
C23 = 2 * (e2.*e3 - e1.*e4);
C31 = 2 * (e3.*e1 - e2.*e4);
C32 = 2 * (e2.*e3 + e1.*e4);
C33 = 1 - 2*e1.^2 - 2*e2.^2;

% Find precession and nutation angle
% Use Body 1-3-1 DCM from Supplementary Materials
kappaAll = acosd(C11);
for i = length(tSpan)
    if (kappaAll(i) > 360)
        kappaAll(i) = kappaAll(i) - 360;
    end
end

phiAll = atan2d(C31,C21);
etaAll = atan2d(-C13,C12);

for i = 1:length(phiAll)
    if (phiAll(i) < 0)
        phiAll(i) = phiAll(i) + 360;
    end
    if (etaAll(i) < 0)
        etaAll(i) = etaAll(i) + 360;
    end
end

% Plot precession vs time
figure;
plot(tSpan, phiAll);
title('Precession Angle Over Time');
xlabel('Time (seconds)');
ylabel('Precession Angle \phi (degrees)');

% Plot nutation vs time
figure;
plot(tSpan, kappaAll);
title('Nutation Angle Over Time');
xlabel('Time (seconds)');
ylabel('Nutation Angle \kappa (degrees)');

figure();
plot(tSpan, e1, tSpan, e2, tSpan, e3, tSpan, e4); grid on;
title("Euler Angle vs. Time");
xlabel("Time [seconds]");
ylabel("Euler Angles [rad]");
legend("E1", "E2", "E3", "E4");

%% ------------ Functions ------------ %%

function xdot = eom(x)

T = 0.1;    % Torque (N-met)
J = 396.2;     
I = 1867;

w1 = x(1);  % Angular Velocity initial conditions
w2 = x(2);
w3 = x(3);
e1 = x(4);  % Euler Parameter initial conditions
e2 = x(5);
e3 = x(6);
e4 = x(7);

% % % Numerical Integration % % %

% Dynamic EOMs
xdot(1) = 0;
xdot(2) = ((I - J) / I) * w1 * w3;
xdot(3) = (T / I) - ((I - J) / I) * w1 * w2;

% Kinematic EOMs
xdot(4) = 0.5 * (w1*e4 - w2*e3 + w3*e2);
xdot(5) = 0.5 * (w1*e3 + w2*e4 - w3*e1);
xdot(6) = 0.5 * (-w1*e2 + w2*e1 + w3*e4);
xdot(7) = -0.5 * (w1*e1 + w2*e2 + w3*e3);

xdot = xdot';

end