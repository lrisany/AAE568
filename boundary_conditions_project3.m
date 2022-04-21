function bcs = boundary_conditions_project3(ya,yb)

I1 = 396.2;
I2 = 1867;
I3 = 1987.8;

J1 = (I3 - I2) / I1;
J2 = (I1 - I3)/I2;
J3 = (I2 - I1)/I3;

L1 = -yb(12)/I1;
L2 = -yb(13)/I2;
L3 = -yb(14)/I3;

lambda1 = yb(12);
lambda2 = yb(13);
lambda3 = yb(14);

lambda = [lambda1;lambda2;lambda3];

w1 = yb(9);
w2 = yb(10);
w3 = yb(11);

G = .5 * [0 -w1 -w2 -w3;w1 0 w3 -w2;w2 -w3 0 w1;w3 w2 -w1 0];

B0 = yb(1);
B1 = yb(2);
B2 = yb(3);
B3 = yb(4);

Beta = [B0;B1;B2;B3];

Betadot = G * Beta;

Gamma1 = yb(5);
Gamma2 = yb(6);
Gamma3 = yb(7);
Gamma4 = yb(8);
gamma = [Gamma1; Gamma2; Gamma3; Gamma4];

w1dot = -J1 * w2 * w3 - lambda1/I1^2;
w2dot = -J2 * w3 * w1 - lambda2/I2^2;
w3dot = -J3 * w1 * w2 - lambda3/I3^2;

omegadot = [w1dot;w2dot;w3dot];

bcs = [ya(1) - 1;     ya(2) - 1 ;  ya(3) - 1; ya(4) - 1;ya(9) - .25 ; ya(10) - .25 ; ya(11) - .25 ;    yb(1);     yb(2);    yb(3);yb(4) ; yb(9); yb(10); yb(11) - .5236;  .5*((L1)^2 + (L2)^2 + (L3)^2) + gamma.' * Betadot + lambda.' * omegadot ];
 
end