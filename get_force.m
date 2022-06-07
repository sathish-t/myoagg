function [force, gbeadfv, gppf, bmmatpt] = get_force(rbead,rmyo,bancm,rdiff,...
    ifor,ipt,bmmat,nhead,fhead,fheadmyp,fone,nsat,kcap,kcap_p,dbead,...
    vmyo0,kod2,r_ring,kwall,kexc,rexc,xmat,rx0,kx)
    
    %{
    calculate the forces on all particles given the velocities
    
    arguments:
        rbead,rmyo - coordinates of actin and myosin beads
        bancm - boolean arrays marking which components are anchored
        rdiff - vector of vectors pointing from a bead to the next along a filament
        ifor - indices of formins in rbead
        ipt - indices of pointed ends in rbead
        bmmat - boolean array showing which actin bead is bound to which myosin
        nhead - number of myosin heads per myosin bead
        fhead,fheadmyp - force per head of myo2 and myp2
        fone - force per filament exerted by myosin if its not saturated with actin
        nsat - threshold number of filaments at which myosin force is saturated
        kcap,kcap_p - spring constants of myo2 and myp2 capture forces
        dbead - distance b/w adjacent actin beads constrained to this value
        vmyo0 - load-free velocity of myo2
        kod2 - actin bending force spring constant
        r_ring - ring radius
        kwall - spring constant for the wall force
        kexc,rexc - myo2-myo2 excluded volume spring constant and range
        xmat - boolean crosslinker matrix with (i,j) = True if ith and jth
            beads are crosslinked
        kx,rx0 - crosslinker spring constant and rest length
        
    returns:
        force - array of forces on all particles
        gbeadfv - array of sum of effective drags from all myosin clusters for each bead
        gppf - effective drag coefficient of each myosin cluster per filament due to the force-velocity
            relationship
        bmmatpt - boolean array with True at (i,j) if ith actin bead is 
            interacting jth myosin and pointed ends are included.
    %}
    
    
    
    %% myosin pulling stall force
    %% ==============================
    % stall force per filament of each myosin cluster, full vector
    sfpf = stall_force_pf (bmmat,fone,nsat,nhead,fhead,fheadmyp,bancm);
    % total stall force amplitude on each bead, full vector
    fsamp = bmmat * sfpf;
    % make it a row vector
    fsamp = fsamp';
    % the contribution from myosin pulling force
    fbead = rdiff .* [fsamp; fsamp; fsamp];
    % sum of rdiff of all beads that interact with each myosin cluster
    rdiffsum = rdiff * bmmat;
        rdiffp = rdiff;
        rdiffp(:,ipt) = rdiffp(:,ipt-1);            % ridffp is the unit tangent vector at all beads including pointed ends
    % the contribution from myosin pulling force
    fmyo = - rdiffsum .* repmat(sfpf',3,1);

    % effective drag coefficient (gamma_pull)of each myosin cluster. NOTE: PER FILAMENT!!!
    gppf = sfpf / vmyo0;
    % the sum of effective drags from all myosin clusters
    gbeadfv = bmmat * gppf;
    % make it a row vector
    gppf = gppf';
    
    %% myosin capture force
    %% ==============================   
    % use bmmatpt which includes pointed end
    bmmatpt = ptcap(rbead, rmyo, ipt, dbead*dbead, bmmat);
    [ib, im] = find(bmmatpt);
    for i = 1:numel(ib)
        this_ib = ib(i);
        this_im = im(i);
        rbm = rmyo(:,this_im) - rbead(:,this_ib);   % r from b to m
        rbm_perp = rbm - sum(rbm .* rdiffp(:,this_ib)) * rdiffp(:,this_ib);
        
        if bancm(this_im)
            this_kcap = kcap;
        else
            this_kcap = kcap_p;
        end
        f_m2b = this_kcap * rbm_perp;
        fbead(:,this_ib) = fbead(:,this_ib) + f_m2b;
        fmyo(:,this_im) = fmyo(:,this_im) - f_m2b;
    end
        
    %% actin filament bending force
    %% ==============================
    % construct fi
    fi = true(1,size(rbead,2));
    fi(ipt) = false;
    fi(ifor) = false;
    fi = find(fi);
    % construct f_i-1
    fimo = true(1,size(rbead,2));
    fimo(ifor) = false;
    fimo(ifor + 1) = false;
    fimo = find(fimo);
    % construct f_i+1
    fipo = true(1,size(rbead,2));
    fipo(ipt) = false;
    fipo(ipt - 1) = false;
    fipo = find(fipo);
    % construct t_i-1, t_i+1 and t_i+2
    tim1 = circshift(rdiff,[0 2]);
    ti = circshift(rdiff,[0 1]);
    tip1 = rdiff;
    tip2 = circshift(rdiff,[0 -1]);
    % bending force
    fb = zeros(size(fbead));
    for i = fimo
        fb(:,i) = kod2 * (tim1(:,i) - (tim1(:,i)' * ti(:,i)) * ti(:,i));
    end
    for i = fi
        fb(:,i) = fb(:,i) + kod2 * (-(1 + ti(:,i)' * tip1(:,i)) * ti(:,i)...
            + (ti(:,i)' * tip1(:,i) + 1) * tip1(:,i));
    end
    for i = fipo
        fb(:,i) = fb(:,i) + kod2 * ((tip1(:,i)' * tip2(:,i)) * tip1(:,i) - tip2(:,i));
    end
    fbead = fbead + fb;
    
    
    %% plasma membrane -- restricting everything inside
    %% ==============================
    % decide which beads / myo go out of the membrane
    radbead = sqrt(sum(rbead(1:2,:) .* rbead(1:2,:)));
    rring = r_ring;
    %rring = sqrt(rring^2 + .05^2)
    outind = find(radbead > rring);
    for i = outind
        fwall = - kwall * rbead(1:2,i) / radbead(i) * (radbead(i) - rring);
        fbead(1:2,i) = fbead(1:2,i) + fwall;
    end

    radmyo = sqrt(sum(rmyo(1:2,:) .* rmyo(1:2,:)));
    outind = find(radmyo > rring);
    for i = outind
        fwall = - kwall * rmyo(1:2,i) / radmyo(i) * (radmyo(i) - rring);
        fmyo(1:2,i) = fmyo(1:2,i) + fwall;
    end
    
    %% myosin excluded volume
    %% ==============================
    % go through every myosin cluster
    for i = 1:size(rmyo,2)
        % difference between the x component of all myosin clusters and
        % the x component of this myosin cluster
        delx = rmyo(1,:) - rmyo(1,i);
        % these are the ones with |delx| < rexc
        ind = abs(delx) < rexc;
        % exclude itself
        ind(i) = false;
        % switch from boolean mode to index mode
        ind = find(ind);
        % difference in y direction 
        dely = rmyo(2,ind) - rmyo(2,i);
        % these are the ones with |delx| < rexc and |dely| < rexc
        ind = ind(abs(dely) < rexc);
        % if no myosin satisfies this criterion, skip to next
        if isempty(ind)
            continue
        end
        % difference in all 3 directions
        delx = delx(ind);
        dely = dely(abs(dely) < rexc);
        delz = rmyo(3,ind) - rmyo(3,i);
        % total distance 
        delr = delx .* delx + dely .* dely + delz .* delz;
        delr = sqrt(delr);
        % these are the ones with delr < rexc
        ind = ind(delr < rexc);
        % if no other myosins within rexc, skip to next 
        if isempty(ind)
            continue
        end
        % delx, dely, delz and delr for those
        delx = delx(delr < rexc);   
        dely = dely(delr < rexc);
        delz = delz(delr < rexc);
        delr = delr(delr < rexc);
        
        % go through every myosin within rexc
        for j = 1:length(ind)
            % magnitude of force
            fmag = kexc * (rexc - delr(j));
            % unit vector of delr
            delrunit = [delx(j); dely(j); delz(j)];
            delrunit = delrunit / norm(delrunit);
            % vector force
            f = - fmag * delrunit;
            % incoporate into fmyo
            fmyo(:,i) = fmyo(:,i) + f;
        end
    end   
        
    %% crosslinker force
    %% ==============================    
    [row, col] = find(xmat);
    for i = 1:numel(row)
        thisrow = row(i);
        thiscol = col(i);
        drv = rbead(:,thisrow) - rbead(:,thiscol);      % vector of dr
        dra = sqrt(sum(drv .* drv));                    % amplitude of dr
        fspa = kx *(rx0 - dra);                         % amplitude of spring force. Positive means compression.
        fspv = fspa * drv / dra;                  % vector of spring force on the bead thisrow.
        fbead(:,thisrow) = fbead(:,thisrow) + fspv;
        fbead(:,thiscol) = fbead(:,thiscol) - fspv;         
    end
    

    % concatenate fbead and fmyo into a 1 by 3*(nbead + nmyo) vector
    force = [fbead, fmyo];
    force = force(:)';
end
