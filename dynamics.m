function dyn = dynamics(x, u)
w1 = x(1);
w2 = x(2);
w3 = x(3);
b1 = x(4);
b2 = x(5);
b3 = x(6);
b4 = x(7);
L1 = u(1);
L2 = u(2);
L3 = u(3);

I1 = 396.2;
I2 = 1867;
I3 = 1987.8;

J1 = (I3 - I2) / I1;
J2 = (I1 - I3) / I2;
J3 = (I2 - I1) / I3;

dyn = [(-J1/I1)*w2*w3 + L1/I1; ...
    (-J2/I2)*w3*w2 + L2/I2; ...
    (-J3/I3)*w1*w2 + L3/I3; ...
    0.5 * (w1*b4 - w2*b3 + w3*b2); ...
    0.5 * (w1*b3 + w2*b4 - w3*b1); ...
    0.5 * (-w1*b2 + w2*b1 + w3*b4); ...
    -0.5 * (w1*b1 + w2*b2 + w3*b3)];
end





