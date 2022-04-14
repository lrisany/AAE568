% disturbance torques
mu_earth = 3.986004418e14;
r_earth = 6371000; 

% TODO: get actual end time AND radius AND theta (earth)
t_endEarth = 86400; % near earth for 1 days (seconds)
t_nearEarth = linspace(1:t_endEarth, 1000); 

r_Earth = linspace(1, 100*r_Earth, 1000); % r is time history of orbit radius from Jup to SC
theta_Earth = linspace(1:t_endEarth, 1000); % angle bt local vertical and SC principle axis

% assume all torques are scalars (will be max value)
% gradity gradient disturbance torque
Tg_Earth = (3*mu_Earth / (2*r_Earth^3)) * abs(I3 - I2) * sin(2*theta_Earth);  % gravity gradiant near earth

% magnetic field disturbance torque
M_earth = 7.96e15; % earth's magnetic moment 
B_earth = M_earth / (r_earth^2);

% sc dipole source: https://commons.erau.edu/cgi/viewcontent.cgi?article=2705&context=space-congress-proceedings
sc_dipole = 4.2 * (0.0001/ (3.335641e-10)); % Ampere*m^2
Tm_earth_man = sc_dipole * B_earth;
Tm_Earth_check = 40e-9; % from https://ntrs.nasa.gov/api/citations/19690020961/downloads/19690020961.pdf

% aerodynamic disturbance torque
Ta_Earth = 0; % assume 0 since little time and far away for all

% solar radiation pressure disturbance torque
% ref: https://www.pveducation.org/pvcdrom/properties-of-sunlight/solar-radiation-in-space
% TODO: get better cmcp distance and area of illuminated surface
sc_rcmcp = 1; % assume center of pressure is 1m away from center of mass (worst case)
sc_As = 1; % area of illuminated surface of SC
q = .8; % reflectance factor for aluminum
psi_earth = 1367; % solar flux density at earth

% find solar flux density as function of AU
% TODO: get vector of radius of SC to sun
R_sun = 695e6; % radius of sun
H_sun = 64e6; % W/m^2 radiant solar intensity at Sun surface
D = linspace(150e9, 150.001e9, 1000); % distance from sun to SC (earth to a earth to a little way away)
psi = ((R_sun^2) / (D^2)) * H_sun;

% TODO: find sun incidence angle over time
i = linspace(0, pi/2, 1000);
c = 3e8; % speed of light
F_sp_earth = (psi / c) * sc_As * (1 + q) * cos(i);

Tsp_ = sc_rcmcp * F_sp_earth;
