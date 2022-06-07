function [rbead,rmyo,xmat,ipt,ifor,bancf,bancm,dbead_first] = update(...
    rbead,rmyo,bmmat,vbead,vmyo,xmat,bancf,bancm,ifor,ipt,kofffor,...
    koffmyo,koffmyp,dt,d_for,d_myo,dbead,dbead_first,rsev,koffx,rxbind,...
    rcap,l_break_sq,r,r_ring,wr,binding_rate_myo2,binding_rate_myp2,...
    binding_rate_for,vpol,binding_rate_x)
    
    %{
    add or remove ring components
    
    arguments:
        rbead,rmyo - coordinates of actin and myosin beads
        bmmat - boolean array showing which actin bead is bound to which myosin
        vbead,vmyo - actin and myosin velocities
        xmat - crosslinker array
        bancf,bancm - boolean arrays marking which components are anchored
        ifor - indices of formins in rbead
        ipt - indices of pointed ends in rbead
        kofffor, koffmyo, koffmyp - respective off rates
        dt - simulation time step
        d_for, d_myo - distances at which formins and myosins are constrained to lie
            from the membrane
        dbead - constraint distance between adjacent actin beads
        dbead_first - vector of distances between formins and next actin bead of each filament
            i.e. between the first and second bead on each filament
        rsev - cofilin mediated severing rate of actin
        koffx - rate of actin crosslinker unbinding
        rxbind - maximum distance between actin beads of different filaments at which
            a crosslinker can bind
        rcap - capture radius of myosin
        l_break_sq - critical length of crosslinker breakage
        r - initial ring radius (irrelevant as the ring does not constrict)
        r_ring - ring radius
        wr - myo and formin bind over a region defined by +- wr along the z axis.
        binding_rate_myo2,binding_rate_myp2,binding_rate_for,binding_rate_x -
            respective binding rates of incoming components
        vpol - rate at which actin filaments grow.
        
    returns:
        rbead,rmyo,xmat,ipt,ifor,bancf,bancm,dbead_first -
            arrays have same meaning as inputs

    %}
    
    % indices of beads
    ibead = 1:length(rbead);
    % make it into vector form
    [xrow,xcol,xval] = find(xmat);
    
    % the default form of bmmat_trial, without any turnover, is equal to
    % bmmat
    bmmat_trial = bmmat;

    %% formin removal
    %% ==============================
    % number of formin dimers to remove
    nforturn = kofffor * dt * length(ifor);
    % make it an integer
    nforturn = floor(nforturn) + (rand < (nforturn - floor(nforturn)));
    % indices of filament to be removed
    ifiloff = randperm(length(ifor),nforturn);
    % sort ifiloff from largest to smallest
    ifiloff = sort(ifiloff,'descend');
    % modify
    for i = ifiloff
        vbead(:,ifor(i):ipt(i)) = [];
        rbead(:,ifor(i):ipt(i)) = [];
        ibead(ifor(i):ipt(i)) = [];
        bancf(ifor(i):ipt(i)) = [];
        bmmat_trial(ifor(i):ipt(i),:) = []; 
    end 
    nbead = length(ibead);
    % update xmat
    % delete entries in xmat that correspond to disappeared beads
    [~,xrow] = ismember(xrow,ibead);
    [~,xcol] = ismember(xcol,ibead);
    xval = xval(and(xrow~=0,xcol~=0));
    xrownew = xrow(and(xrow~=0,xcol~=0));
    xcol = xcol(and(xrow~=0,xcol~=0));
    xrow = xrownew;
    % update ifor and ipt
    ifor(~ismember(ifor,ibead)) = [];
    ipt(~ismember(ipt,ibead)) = [];
    [~,ifor] = ismember(ifor,ibead);
    [~,ipt] = ismember(ipt,ibead);
    dbead_first(ifiloff) = [];
        
    %% formin binding and nucleation
    %% ==============================
    nfor_in = binding_rate_for * 2 * pi * r_ring * dt;                     % number of formin dimers to nucleate
    nfor_in = floor(nfor_in) + (rand < (nfor_in - floor(nfor_in)));     % make it an integer
    for i = 1:nfor_in
        theta = 2 * pi * rand;
        
        [x, y] = pol2cart(theta,r_ring-d_for);       % obtain the position of the new formin
        z = wr * (rand - .5);
        r_new_for = [x;y;z];
        
        % now pick a random direction in 3D, see
        % http://mathworld.wolfram.com/SpherePointPicking.html, Eq. 16
        x = randn;
        y = randn;
        z = randn;
        vrand = [x;y;z];
        vrand = vrand / norm(vrand);
        if vrand(1:2)' * vrand(1:2) > r_ring*r_ring
            vrand = - vrand;     % cannot nucleate into the membrane       
        end
        r_new_act = r_new_for + 0.1 * dbead * vrand;    % position of the new actin bead
        rbead = [rbead, r_new_for, r_new_act];
        ipt = [ipt, ipt(end) + 2];
        ifor = [ifor, ipt(end) - 1];
        nbead = nbead + 2;
        bancf = [bancf, true, false];
        bmmat_trial = [bmmat_trial; false(2,size(bmmat_trial,2))];
        dbead_first = [dbead_first, 0.1 * dbead];
    end
        
    %% cofilin severing
    %% ==============================
    % total length of actin filament
    ltot = dbead * (nbead - length(ifor));
    % number of severing events
    nsev = ltot * rsev * dt;
    % make it an integer
    nsev = floor(nsev) + (rand < nsev - floor(nsev));
    % cannot sever at formin beads or the bead after (leaving a bare formin)
    legalsite = 1:nbead;
    legalsite(or(ismember(legalsite,ifor),ismember(legalsite-1,ifor))) = [];
    % sites of severing
    isev = legalsite(randperm(length(legalsite),nsev));
    % sort isev from largest to smallest
    isev = sort(isev,'descend');

    % modify
    ibead = 1:size(rbead,2);
    bfor = ismember(ibead,ifor); % boolean vector of formins
    bpt = ismember(ibead,ipt);
    for i = isev
        % the pointed end of the severed filament
        iptsev = find(bpt(i:end), 1) + i - 1; % note: ipt must be already sorted from smallest to largest
        vbead(:,i:iptsev) = [];
        rbead(:,i:iptsev) = [];
        bmmat_trial(i:iptsev,:) = [];   
        ibead(i:iptsev) = [];
        bancf(i:iptsev) = [];
        bfor(i:iptsev) = [];
        bpt(i:iptsev) = [];
        bpt(i-1) = true;
    end
    ifor = find(bfor);
    ipt = find(bpt);
    
    nbead = size(rbead,2); 
    
    % update xmat
    % delete entries in xmat that correspond to disappeared beads
    [~,xrow] = ismember(xrow,ibead);
    [~,xcol] = ismember(xcol,ibead);
    xval = xval(and(xrow~=0,xcol~=0));
    xrownew = xrow(and(xrow~=0,xcol~=0));
    xcol = xcol(and(xrow~=0,xcol~=0));
    xrow = xrownew;
         
        
    %% myosin turnover
    %% ==============================
    % number of myosin clusters to remove
    nmyoturn = koffmyo * sum(bancm) * dt;
    nmypturn = koffmyp * sum(~bancm) * dt;
    % make it an integer
    nmyoturn = floor(nmyoturn) + (rand < nmyoturn - floor(nmyoturn));
    nmypturn = floor(nmypturn) + (rand < nmypturn - floor(nmypturn));
    % indices of myosin clusters to remove
    temp = find(bancm);
    imyooff = randperm(sum(bancm), nmyoturn);
    imyooff = temp(imyooff);
    temp = find(~bancm);
    imypoff = randperm(sum(~bancm), nmypturn);
    imypoff = temp(imypoff);        
    imyooff = [imyooff, imypoff];
    % modify
    vmyo(:,imyooff) = [];
    rmyo(:,imyooff) = [];
    bancm(imyooff) = [];
    bmmat_trial(:,imyooff) = []; 
    
    nmyo2_in = binding_rate_myo2 * 2 * pi * r_ring * dt;                % number of myosin clusters to come in
    nmyo2_in = floor(nmyo2_in) + (rand < nmyo2_in - floor(nmyo2_in));   % make it an integer
    for i = 1:nmyo2_in
        theta = 2 * pi * rand;
        if (rand() > abs(sin(theta)))
            theta = 2 * pi * rand;
        end
        
        [x_in, y_in] = pol2cart(theta,r_ring-d_myo);
        z_in = wr * (rand - .5);
        rmyo = [rmyo,[x_in; y_in; z_in]];
        
        bancm = [bancm, true];
        bmmat_trial = [bmmat_trial, false(size(bmmat_trial,1),1)]; 
    end
    
    nmyp2_in = binding_rate_myp2 * 2 * pi * r_ring * dt;                           % number of myosin clusters to come in
    nmyp2_in = floor(nmyp2_in) + (rand < nmyp2_in - floor(nmyp2_in));   % make it an integer
    for i = 1:nmyp2_in
        ind = 1:size(rbead,2);
        ind(ifor) = [];                     % indices for actin beads (non-formin)
        ind = ind(randi(numel(ind)));  % pick one at random
        for itrial = 1:1e3
            x = randn;      % pick a random direction in 3D
            y = randn;
            z = randn;
            vrand = [x;y;z];
            vrand = vrand / norm(vrand);
            % position of new myp2
            r_new = rbead(:,ind) + rcap * vrand;
            if r_new(1:2)' * r_new(1:2) < r_ring*r_ring   % stop trying if some r_new that's within the membrane is found
                break
            end
        end
        rmyo = [rmyo,r_new];
        bancm = [bancm, false];
        bmmat_trial = [bmmat_trial, false(size(bmmat_trial,1),1)];
    end
    
    
    %% crosslinker breaking due to over extension
    %% ==============================
    cross_dx = rbead(:,xrow) - rbead(:,xcol);
    cross_dx2 = sum(cross_dx .* cross_dx);
    ind_del = cross_dx2 > l_break_sq;
    xrow(ind_del) = [];
    xcol(ind_del) = [];
    xval(ind_del) = [];
    xmat = sparse(xrow,xcol,xval,nbead,nbead);
    
    
    %% crosslinker turnover in itself
    %% ==============================
    % number of xlinkers to go off
    nxturn = koffx * dt * .5 * full(sum(sum(xmat)));
    % make it an integer
    nxturn = floor(nxturn) + (rand < nxturn - floor(nxturn));
    % random indices to go off
    ixoff = randperm(sum(sum(xmat))/2,nxturn);
    % delete these xlinkers
    [xrow,xcol,xval] = find(triu(xmat));
    xval(ixoff) = false;
    xmat = sparse(xrow,xcol,xval,nbead,nbead);
    xmat = or(xmat,xmat');
    % number of crosslinkers after breaking
    nafter = sum(xmat(:)) * .5;
    % add a number of xlinkers such that the total is constant
    xmat = full(xmat);
    nadd = binding_rate_x * 2 * pi * r_ring * dt;
    nadd = floor(nadd) + (rand < nadd - floor(nadd));
    for iadd = 1:nadd
        for itry = 1:1000
            if itry == 1000
                error('too many trials')
            end
            % pick a random bead
            thisxrow = randi(nbead);
            % find a neighbor of this bead
            xnbr = abs(rbead(1,:) - rbead(1,thisxrow)) < rxbind;
            xnbr(thisxrow) = 0;
            candidate = and(xnbr, abs(rbead(2,:) - rbead(2,thisxrow)) < rxbind);
            % this bead and candidate cannot be already linked
            candidate(xmat(thisxrow,:)) = false;
            
            candidate = find(candidate);
            dx = rbead(1,candidate) - rbead(1,thisxrow);
            dy = rbead(2,candidate) - rbead(2,thisxrow);
            dz = rbead(3,candidate) - rbead(3,thisxrow);
            candidate = candidate(dx.*dx+dy.*dy+dz.*dz < rxbind * rxbind);

            % if no neighbor can be found
            if isempty(candidate)
                continue
            else
                thisxcol = randi(length(candidate));
                thisxcol = candidate(thisxcol);
                xmat(thisxrow,thisxcol) = true;
                xmat(thisxcol,thisxrow) = true;
                break
            end
        end
    end
        
        
    %% actin polymerization (add bead)
    %% ==============================
    [xrow,xcol,xval] = find(xmat);
    poly_ind = find(dbead_first > 1.1 * dbead);
    for ipoly = numel(poly_ind) : -1 : 1            % reverse order, convenient for insersion operations
        i = poly_ind(ipoly);
        dbead_first(i) = dbead_first(i) - dbead;
        r1 = rbead(:,ifor(i));
        r2 = rbead(:,ifor(i)+1);
        dr = r2 - r1;
        r_new = r2 - dbead * dr / norm(dr);     % position of the new bead to add
        rbead = [rbead(:,1:ifor(i)), r_new, rbead(:,ifor(i)+1:end)];
        vbead = [vbead(:,1:ifor(i)), [0;0;0], vbead(:,ifor(i)+1:end)];
        bancf = [bancf(1:ifor(i)), false, bancf(ifor(i)+1:end)];
        bmmat_trial = [bmmat_trial(1:ifor(i),:); false(1,size(bmmat_trial,2)); bmmat_trial(ifor(i)+1:end,:)];
        nbead = nbead + 1;
        ifor = ifor + 1;
        ifor(1:i) = ifor(1:i) - 1;              % these two lines are equivalent to ifor(i+1:end) += 1, but avoid the endpoint problem
        ipt(i:end) = ipt(i:end) + 1;
        xrow = xrow + (xrow > ifor(i));         % to update xmat, it is sufficient to update xrow and xcol
        xcol = xcol + (xcol > ifor(i));
    end
    xmat = sparse(xrow,xcol,xval,nbead,nbead);
    dbead_first = dbead_first + vpol * dt;
    
end
