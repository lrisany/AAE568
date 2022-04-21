function dydt = BVP_ode_project3(t,y)

I1 = 396.2;
I2 = 1867;
I3 = 1987.8;

J1 = (I3 - I2) / I1;
J2 = (I1 - I3)/I2;
J3 = (I2 - I1)/I3;

B0 = y(1,:);
B1 = y(2,:);
B2 = y(3,:);
B3 = y(4,:);

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



G = .5 * [0 -w1 -w2 -w3;w1 0 w3 -w2;w2 -w3 0 w1;w3 w2 -w1 0];

gamma = [Gamma1; Gamma2; Gamma3; Gamma4];

omega = [w1;w2;w3];

Beta = [B0; B1; B2; B3];

Betadot = G * Beta;
Gammadot = G * gamma;

lambdadot = [0, J3*lambda3, J2*lambda2; J3 *lambda3, 0, J1*lambda1;J2*lambda2, J1*lambda1,0] * omega + .5 * [B1, -B0, -B3, B2;B2, B3,-B0, -B1;B3,-B2,B1,-B0] * gamma;
w1dot = -J1 * w2 * w3 - lambda1/I1^2;
w2dot = -J2 * w3 * w1 - lambda2/I2^2;
w3dot = -J3 * w1 * w2 - lambda3/I3^2;

omegadot = [w1dot;w2dot;w3dot];


dydt = y(15) * [Betadot; Gammadot ; omegadot;lambdadot ; 0];

end