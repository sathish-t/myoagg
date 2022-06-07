%{ 
   set up initial ring configuration and some more parameters
   for the simulation 
%}


% Geometry
% radius of the intact cell (3.7 / 2 in septation paper)
% ST: increasing this to 2x.
r = 3.7;      
lr = r * 2 * pi;    
lf = 1.3; % mean length of actin filament in the ring
% Initial anchored region(s)
% la is a 2 by n matrix, each column specifies the start and end of an
% anchored region
la = lr * [0; 1];

% binding radius of ain1
rxbind = 0.05;
% width of the ring
wr = 0.2;
% depths of formin and Myo2 clusters under the membrane
d_for = .03;        % 0.03 micron from sup res paper Fig. 2F
d_myo = .08;        % 0.08 micron from sup res paper Fig. 2D
% distance between beads on an actin filament
dbead = 0.1;
dbead_sq = dbead ^2;
% capture radius of myosin clusters (unfortunately the name rmyo is
% used by another variable)
rcap = 0.08; % if this is modified, remember to edit ptcap.m
% minimum distance between myosin clusters (excluded volume)
rexc = 0;

% density of components
nsp = 2.5;
rhom = 7.5 * nsp;
rhof = 15;
rhoain = 25;
    
for itrial = 1:1e5
    
    %% Initial positions of all molecules
    [r, rbead, rmyo, ipt, ifor, bancf, bancm, xmat] = initial_circle_rand (rhom, rhof, rhoain, lr, ...
        la, wr, d_for, d_myo, dbead, rexc, rxbind);
    fil_tag = 1:numel(ipt);
    r_ring = r; 
    confbad = false;
    % test if any two myosin clusters are within excluded radius
    for i = 1:size(rmyo,2)
        xtemp = abs(rmyo(1,i)-rmyo(1,:)) < rexc;    % myosins that are within rexc of the current myosin in the x direction
        ytemp = abs(rmyo(2,i)-rmyo(2,:)) < rexc;
        ztemp = abs(rmyo(3,i)-rmyo(3,:)) < rexc;
        xytemp = and(xtemp, ytemp);
        xyztemp = and(xytemp,ztemp);
        xyztemp(i) = false;                         % now xyztemp are myosins within rexc of the current myosin in x,y and z directions
        if any(xyztemp)     % if there is any such candidate
            for j = find(xyztemp)
                rmyo_diff = rmyo(:,j) - rmyo(:,i);
                rmyo_diff_sq = rmyo_diff' * rmyo_diff;
                if rmyo_diff_sq < rexc * rexc
                    confbad = true;
                    break
                end
            end
        end
        if confbad
            break
        end
    end
    if ~confbad
        break
    end
end

%% if a good configuration is found, output
if ~confbad
    vbead = zeros(size(rbead));
    vmyo = zeros(size(rmyo));
else
    error('all trials gave bad configurations')
end

% number of myosin clusters
nm = size(rmyo,2);