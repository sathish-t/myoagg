%{ 
   set up some parameters for the simulation 
%}

% constraints are enforced over a timescale tau_o_dt * dt
% where dt is the simulation timestep
tau_o_dt = 10;            

% Forces
% Static force of a myosin head on a filament (pN) (1 in exp)
fhead = 2.5;
fheadmyp = 2.5;
% number of heads per myosin cluster
nhead = 40 / nsp;
% maximal stall force per myosin cluster per filament
fone = 4;
% Maximum number of filaments that a myosin cluster interacts with, without overloading (10 in dev cell)
nsat = fhead*nhead/fone;

% load-free velocity
vmyo0 = 0.24;
% spring constant for myosin - actin capture force
kcap = 400%40 / nsp or 100;
kcap_p = kcap/4;

% Myosin cluster drag coefficient from membrane
gm = 1300 / nsp * 3 *6;
% Formin drag coefficient from membrane
gf = 1900;
% membrane elastic constant
kwall = 20;
% breaking length of xlinkers (.05 in dev cell)
l_break = .05;
l_break_sq = l_break^2;    
% peeling force
kexc = 10000 / nsp;    % myo-myo excluded volume spring constant
rexc = 45/1000; % myo-myo xvol range
% artificial drag between adjacent actin beads on the same filament
ga = 4;

%% constraints
% septation rate (absolute value; defined as dr/dt; .07 micron / minute in Pelham and Chang)
% set to zero here
vs_abs = .07/60 * 0
% crosslinker spring constant (25 in Bala. for a-actinin)
kx = 25;
% corsslinker rest length
rx0 = 0.03;

% water viscosity (in pN um-2 / s)
mu = 0.6;
% solution drag of actin beads, see Broersma 1960
half_length = dbead / 2;
dact = 0.010;            % diameter of actin filament
b = dact / 2;
sigma = log(2*half_length / b);
gamma = 0.35 - 4 * (1/sigma - 0.43)^2;
gb = 4 * pi * mu * dbead / (sigma - gamma);
% solution drag of myosin clusters 
gmsol = 6 * pi * mu * rcap;
% bending stiffness of actin filament
kappa = 0.041;
% kappa over dbead^2
kod2 = kappa / (dbead^2);
% maximal bending force on a bead
maxfb = 1e6; % 0.14 at 3kT

% turnover
% Actin severing rate per filament length by cofilin (46% actin left
% after 60s of constriction. cf. 1.8 / 60 um-1 s-1 in dev cell)
rsev = 1.8/60;
% take-out rate of formins (73% actin left after 60s of constriction with Jasp. cf. 0.023 in dev cell.)
kofffor = 0.023; 
% take-out rate of myosins (79% left after 60s of constriction in Bala. cf. 0.026 in dev cell.)
koffmyo = 0.026;
koffmyp = 0*koffmyo;
% take-out rate of a-actinin (3.3 in Dev Cell.)
koffx = 3.3;

% binding rates of myosins, formin and xlinkers
binding_rate_myo2 = koffmyo * rhom; %(0.2 in dev cell)
binding_rate_myp2 = 0* koffmyp * rhom * 2/3;
binding_rate_for = kofffor * rhof; % 0.35 in dev cell
binding_rate_x = koffx * rhoain;

% actin growth rate
vpol = 0.07; % (0.07 in dev cell)
dbead_first = repmat(dbead,1,numel(ifor));

% initialize velocity and constraint force arrays
vbead = 0 * rbead;
vmyo = 0 * rmyo;
fc = zeros(3*(size(rbead,2)+size(rmyo,2)),1);
poly_time = false;