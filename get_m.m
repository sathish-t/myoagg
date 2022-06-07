function w = get_m (ga,ifor,ipt,rbead,rmyo,bancf,bancm,bmmat,rdiff,gbeadfv,gppf,gb,gf,gm,gmsol,bmmatpt)

    %{ 
    calculates the drag coefficients of all components as a 2D matrix.
    matrix contributions come from: force velocity relationship of myosin pulling forces,
    solution drag on actin and myosin beads, membrane drag on myo2 and formin beads, artificial
    drag between adjacent actin beads and between actin bead and bound myosins.
    
    NOTE: for the artifical drags, only the diagonal elements are calculated here. off-diagonal
    elements are calculated separately in get_a.
    
    First w of each individual particle is calculated, then concatenated in
    block-diagonal form.
    
    arguments:
        ga - artificial drag between actin beads on the same filament, and between
            an actin bead and its bound myosin
        ifor - indices of beads corresponding to formins
        ipt - indices of beads corresponding to actin pointed ends
        rbead, rmyo - coordinates of actin and myosin beads
        bancf,bancm - boolean vector of anchored actin beads and myosins
        bmmat - boolean matrix where i,j = True if ith actin is bound to jth myosin, excluding pointed ends.
        rdiff - vector of vectors pointing from one act bead to the next along a filament
        gbeadfv - 
        gppf - vector of drag coefficients associated with myosin pulling forces
        gb - drag coefficient of actin beads due to interaction with the solution
        gf - Formin drag coefficient from membrane
        gm - Myosin cluster drag coefficient from membrane
        gmsol - drag coefficient of actin beads due to interaction with the solution
        bmmatpt - same as bmmat but including pointed ends
        
    returns:
        w: the required matrix.

    %}
    
    nbead = size(rbead,2);
    % preallocate
    % total number of particles times 3
    tottimes3 = 3*(size(rbead,2)+size(rmyo,2));
    % global values of row, col and v
    glbrow = [];
    glbcol = glbrow;
    glbv = glbrow;

    %% for beads
    %% ==============================
    totdrag_vec = ga * sum(bmmatpt,2);
    for i = 1:size(rbead,2)
        % cytosolic drag
        m = zeros(3,3);
        m(1,1) = gb;
        m(2,2) = gb;
        m(3,3) = gb;
        % membrane drag
        if bancf(i) % if the ith bead is an anchored formin
            m = m + gf * eye(3);
        end
        % effective drag from f-v relation
        nn = rdiff(:,i) * rdiff(:,i)';
        m = m + (gbeadfv(i)-totdrag_vec(i)) * nn;
        %% artificial drag between adjacent actin beads on each filament
        totdrag = 2 - 2 * any(ifor == i) - any(ifor+1 == i) - any(ipt == i);
        totdrag = ga * totdrag;
        m(1,1) = m(1,1) + totdrag;
        m(2,2) = m(2,2) + totdrag;
        m(3,3) = m(3,3) + totdrag;
        %% artificial drag between myosin and actin beads
        totdrag = totdrag_vec(i);
        m = m + totdrag * eye(3);
        %% 
        % get w
        w = m;
        % store
        % find row, col and v for this w itself (local)
        [locrow, loccol, locv] = find(w);
        % change it to the global form and store
        glbrow = [glbrow; locrow + 3*(i-1)];
        glbcol = [glbcol; loccol + 3*(i-1)];
        glbv = [glbv; locv];
    end

    %% for myosin
    %% ==============================
    % calculate nn for each bead
    nn = zeros(3,3,size(rbead,2));
    for i = 1:size(rbead,2)
        nn(:,:,i) = rdiff(:,i) * rdiff(:,i)';        
    end
    
    for i = 1:size(rmyo,2)
        % solution drag 
        m = gmsol * eye(3);
        % membrane drag
        if bancm(i) % if the ith myosin is anchored
            m = gm * eye(3);
        end
        % sum of nn for all beads interacting with the ith myosin
        nni = sum(nn(:,:,bmmat(:,i)),3);
        % effective drag from f-v relation
        m = m + (gppf(i)-ga) * nni;
        %% artificial drag between myosin and actin beads
        totdrag = ga * sum(bmmatpt(:,i));
        m = m + totdrag * eye(3);
        %% 
        % get w
        w = m;
        % store
        % find row, col and v for this w itself (local)
        [locrow, loccol, locv] = find(w);
        % change it to the global form and store
        glbrow = [glbrow; locrow + 3*(i-1+nbead)];
        glbcol = [glbcol; loccol + 3*(i-1+nbead)];
        glbv = [glbv; locv];
    end
    w = sparse(glbrow,glbcol,glbv,tottimes3,tottimes3);
end