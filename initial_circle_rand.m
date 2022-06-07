function [r, rbead, rmyo, ipt, ifor, bancf, bancm, xmat] = initial_circle_rand(rhom, rhof,...
            rhoain, lr, la, wr, d_for, d_myo, dbead, rexc, rxbind)
    
    %{ 
    set up an initial ring configuration
    
    arguments:
        rhom, rhof, rhoain - densities of myo2, formin, alpha-actinin
        lr - total length of the ring
        la - array representing the section of the ring that is anchored.
            (irrelevant for this simulation as we set it = [0 1] * lr i.e. the entire
            ring is anchored)
        wr - myo and formin bind over a region defined by +- wr along the z axis
        d_for, d_myo - distances at which formins and myosins are constrained to lie
            from the membrane
        dbead - constraint distance between adjacent actin beads
        rexc - range of myo2-myo2 excluded volume
        rxbind - maximum distance between actin beads of different filaments at which
            a crosslinker can bind
        
    returns:
        r - ring radius
        rbead, rmyo - coordinates of actin and myosin beads
        ipt - indices of beads corresponding to actin pointed ends
        ifor - indices of beads corresponding to formins
        bancf,bancm - boolean arrays marking which components are anchored
        xmat - crosslinker array

    %}
    
    %% set the ring shape
    %% ==============================
    % radius of the ring
    r = lr / (2*pi);  
    % square of r
    rsq = r*r;
    % gamma as a tuning angle
    gamma = pi/12;
    % rotation matrix for small gamma and -gamma
    rot1 = [cos(gamma), -sin(gamma); sin(gamma), cos(gamma)];
    rot2 = [cos(gamma), sin(gamma); -sin(gamma), cos(gamma)];
    % number of formin dimers
    nf = round(rhof * lr);
    %% form the base shape of the ring
    % number of points on the spine
    nsp = 1000 * nf;
    % r and theta coordinates of the spine, first as a perfect circle
    rsp = r * ones(1,nsp); % sp means spine
    thsp = linspace(0, 2*pi, nsp +1);
    thsp(end) = [];
    % thanc is a 2 by n matrix, each column specifies the start and end
    % points of an anchored region
    thanc = la / lr * 2 * pi;
    % elements in thsp that corresponds to detached region
    for i = 1:size(la,2)
        if size(la,2) == 1
            ancreg = and(thsp > thanc(1), thsp < thanc(2));
            ind = find(~ancreg);
        elseif i == 1
            ind = and(thsp > thanc(2,1), thsp < thanc(1,2));
            ind = find(ind);
        elseif i == 2
            ind = thsp > thanc(2,2);
            ind = find(ind);
        end
        % find the head of the detached region
        if min(ind) == 1 && max(ind) == nsp
            for j = 2:length(ind)
                if ind(j) ~= ind(j-1) + 1
                    head = j;
                    break
                end
            end
        else
            head = 1;
        end
    end
    % transform into x and y
    [xsp, ysp] = pol2cart(thsp, rsp);
    % tangents at each point
    tgx = xsp - circshift(xsp,[0 1]);
    tgy = ysp - circshift(ysp,[0 1]);
    % normalize
    for i = 1:length(tgx)
        seglength = norm([tgx(i), tgy(i)]);
        tgx(i) = tgx(i) / seglength;
        tgy(i) = tgy(i) / seglength;
    end
    
    
    %% positions of myosin clusters
    %% ==============================
    
    %% output rmyo
    % number of myosin clusters
    nm = round(rhom * lr);
    % x and y coordinates
    % x and y directions: follow spine
    dsp = lr / nsp;                 % distance between adjacent spine points
    nd_min = round(1.1*rexc / dsp);     % minimal separation in terms of spine points
    free_l = nsp - nd_min * nm;     % free length in which myosins can be placed
    ind_free = randperm(free_l,nm); % indices within free_l
    ind_free = sort(ind_free);
    ind = ind_free + nd_min * (0:nm-1);
    xm = xsp(ind);
    ym = ysp(ind);
    % z direction
    zm = -.5*wr + wr * rand(1,nm);
    rmyo = [xm; ym; zm];                                        
    r = lr / 2 / pi;
    rmyo(1:2,:) = rmyo(1:2,:) * (1 - d_myo /r); 
    % cylindrical coordinate
    [thm, ~] = cart2pol(xm,ym);
    thm(thm<0) = thm(thm<0) + 2*pi;
    
    
    %% output rbead
    %% ==============================
    % number of formin dimers                     
    nf = round(rhof * lr);
    % positions of formin dimers
    % z direction
    zf = -.5*wr + wr * rand(1,nf);
    % x and y directions: follow spine
    ind = randperm(nsp,nf);
    xf = xsp(ind) * (1 - d_for /r);
    yf = ysp(ind) * (1 - d_for /r);
    % rf in CYLINDRICAL coordinates
    [thf, ~] = cart2pol(xf,yf);
    thf(thf<0) = thf(thf<0) + 2*pi;
    % probability for a filament to have length lf
    % consider filaments that are dbead to 100dbead long
    lfvec = 1:100;
    lfvec = lfvec * dbead; % 
    % probability 
    aa = (.023 + lfvec * 0.03) .* exp(-lfvec .* (2*.023 + lfvec * .03) / (2 * .07));
    aa = aa / sum(aa); % normalize the total to 1
    % cumulative probability
    acum = cumsum(aa);
    % preallocate
    rbead = [];
    ifor = [];
    ipt = [];
    % go through each filament
    for i = 1:nf
        % put formin
        rbead = [rbead, [xf(i); yf(i); zf(i)]]; % use cylindrical coordinate for now
        % set ifor
        ifor = [ifor, size(rbead,2)];
        % number of actin beads
        nb = find(acum > rand, 1);
        % set ipt
        ipt = [ipt, ifor(end) + nb];
        % polarity
        pol = randi(2)*2 - 3; % a random number, either 1 or -1
        % go through each actin bead
        for j = 1:nb
            % find the spine point closest to last bead
            [~, isp] = min(abs(xsp - rbead(1,end)) + abs(ysp - rbead(2,end)));
            % coordinate for this bead
            xbead = rbead(1,end) + pol * dbead * tgx(isp);
            ybead = rbead(2,end) + pol * dbead * tgy(isp);
            % if this bead sticks out of the ring
            if xbead*xbead + ybead*ybead > rsq
                % try to tune in one direction
                tg = rot1 * [tgx(isp); tgy(isp)];
                xbead = rbead(1,end) + pol * dbead * tg(1);
                ybead = rbead(2,end) + pol * dbead * tg(2);
                if xbead*xbead + ybead*ybead > rsq
                    % try to tune in the other direction
                    tg = rot2 * [tgx(isp); tgy(isp)];
                    xbead = rbead(1,end) + pol * dbead * tg(1);
                    ybead = rbead(2,end) + pol * dbead * tg(2);
                    % if neither works, report error
                    if xbead*xbead + ybead*ybead > rsq
                        %error('nothing works!')
                    end
                end
            end
            % put coordinates into rbead
            rbead = [rbead, [xbead; ybead; zf(i)]];
        end
    end
    nbead = size(rbead,2);
    
    %% output bancf (boolean vector of anchored formins) and bancm (anchored myosin clusters)
    %% ==============================
    bf = zeros(1,nf);
    bancm = zeros(size(thm));
    % anchored region(s)
    % thanc is a 2 by n matrix, each column specifies the start and end
    % points of an anchored region
    thanc = la / lr * 2 * pi;
    
    for i = 1:size(thanc,2)
        bf = bf + and(thf > thanc(1,i),thf < thanc(2,i));
    end
    bf = logical(bf);
    % change bf (formins only) to bancf (all beads)
    bancf = false(1,nbead);
    bancf(ifor) = bf;
    % myosin
    for i = 1:size(thanc,2)
        bancm = bancm + and(thm > thanc(1,i),thm < thanc(2,i));
    end
    bancm = logical(bancm);
    
    %% output xmat (nbead by nbead square boolean matrix showing crosslinking)
    %% ==============================
    % i and j cannot belong to the same filament, this is
    % ensured by setting bead separation larger than maximal
    % length of ain1. check this.
    if rxbind > dbead
        error ('rxbind > dbead')
    end
    % preallocate xmat
    xmat = false(nbead);
    % maximal theta between two beads that can be bound by ain1
    thmax = rxbind / r;
    % theta of rbead
    [thb, ~] = cart2pol(rbead(1,:), rbead(2,:));
    % find all possible pairs within connection length of ain1
    for i = 1:nbead-1
        for j = i+1:nbead
            if norm(rbead(:,i) - rbead(:,j)) < rxbind && ~ismember(i,ifor) && ~ismember(j,ifor)
                xmat(i,j) = true;
            end
        end
    end
    % total number of possible binding sites
    nsite = sum(sum(xmat));
    % total number of ain1
    nain = round(rhoain * lr);
    % there must be enough sites
    if nsite < nain
        error('nsite < nain')
    end
    % indices of xmat that has entry 1
    xind = 1:(nbead*nbead);
    xind = xind(xmat);
    % randomly select nain entries from xind
    xind = xind(randperm(length(xind),nain));
    % clear xmat and put these entries in
    xmat = false(nbead);
    xmat(xind) = true;
    % make it symmetric
    xmat = or(xmat,xmat');
    % make it a sparse matrix
    xmat = sparse(xmat);
end