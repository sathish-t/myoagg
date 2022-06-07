function amat = get_a(rbead, rmyo, bmmat, gppf, rdiff, ga, ifor, ipt, bmmatpt)

    %{
    get the sparse A drag coefficient matrix, which corresponds to the off-diagonal
    elements of the artificial drag between actin beads and myosins to which
    they are bound, and that between adjacent actin beads.
    
    arguments:
        rbead, rmyo - coordinates of actin and myosin beads
        bmmat - boolean matrix where i,j = True if ith actin is bound to jth myosin, excluding pointed ends.
        gppf - vector of drag coefficients associated with myosin pulling forces
        rdiff - vector of vectors pointing from one act bead to the next along a filament
        ga - artificial drag between actin beads on the same filament, and between
            an actin bead and its bound myosin
        ifor - indices of beads corresponding to formins
        ipt - indices of beads corresponding to actin pointed ends
        bmmatpt - same as bmmat but including pointed ends
        
    returns:
        amat: the required matrix.

    %}
    
    
    %% preallocation
    % Estimate: there are nmyo myosin clusters, each interacts with ~20
    % filaments, so ~ 20*nmyo pairs. Each pair gives 18 non-zero entries in
    % the A matrix, so 360 * nmyo non-zero elements. There are nbead beads,
    % each gives 18 non-zero entries due to artificial drag, so 18 * nbead
    % non-zero elements. Use a larger number for safety.
    arow = zeros(1,720*size(rmyo,2));
    acol = arow;
    aval = arow;
    ind = 1;
    indvec = 1:9;
    %% f-v relation 
    % bead-list and myosin-list of interacting bead-myosin pairs
    [blist, mlist] = find(bmmat);
    for i = 1:length(blist)
        ib = blist(i);
        im = mlist(i) + size(rbead,2);
        g = gppf(mlist(i)) - ga; % gppf means gamma_pull per filament. It is a vector with nmyo elements.
        rd = rdiff(:,ib); % it is a column-vector
        % get the local gamma (3 X 3 matrix). It will be put into the sparse matrix A. 
        % for the equation, see notes "effective drag force"
        gmat = g * (rd * rd'); % note that gmat is symmetric
        % put gmat into the sparse matrix A 
        arow(indvec) = [1 2 3 1 2 3 1 2 3] + 3 * (ib-1);
        acol(indvec) = [1 1 1 2 2 2 3 3 3] + 3 * (im-1);
        aval(indvec) = gmat(:);
        indvec = indvec + 9;
        arow(indvec) = [1 2 3 1 2 3 1 2 3] + 3 * (im-1);
        acol(indvec) = [1 1 1 2 2 2 3 3 3] + 3 * (ib-1);
        aval(indvec) = gmat(:);
        indvec = indvec + 9;
        % update ind
        ind = ind + 18;
    end

    %% artificial drag between myosin and actin
    garow = zeros(1,720*size(rmyo,2));
    gacol = garow;
    gaval = garow;
    ind2 = 1;
    gmat = ga * eye(3);
    [grow, gcol, gval] = find(gmat);
    [bptlist, mptlist] = find(bmmatpt);
    for i = 1:length(bptlist)
        ib = bptlist(i);
        im = mptlist(i) + size(rbead,2);
        n = length(grow);
        % put gmat into the sparse matrix ga
        garow(ind2:ind2+n-1) = grow' + 3 * (ib-1);
        gacol(ind2:ind2+n-1) = gcol' + 3 * (im-1);
        gaval(ind2:ind2+n-1) = gval';
        garow(ind2+n:ind2+2*n-1) = grow' + 3 * (im-1);
        gacol(ind2+n:ind2+2*n-1) = gcol' + 3 * (ib-1);
        gaval(ind2+n:ind2+2*n-1) = gval';
        % update ind
        ind2 = ind2 + 2*n;
    end
    %% artifical drag between adjacent actin beads
    for i = 1:size(rbead,2)
        % if this bead is not a formin or the next bead after a formin
        if ~ or(any(ifor == i), any(ifor+1 == i))
            % artificial drag between this bead and the previous bead
            arow(ind:ind+2) = (1:3) + 3*(i-1);
            acol(ind:ind+2) = (1:3) + 3*(i-2);
            aval(ind:ind+2) = repmat(ga,1,3);
            ind = ind + 6;
        end
        % if this bead is not a pointed end or a formin
        if ~or(any(ipt == i), any(ifor == i))
            % artifical drag between this bead and the next bead
            arow(ind:ind+2) = (1:3) + 3*(i-1);
            acol(ind:ind+2) = (1:3) + 3*(i);
            aval(ind:ind+2) = repmat(ga,1,3);
            ind = ind + 6;
        end
    end
    %% gather all entries together
    aval(arow==0) = [];
    arow(arow==0) = [];
    acol(acol==0) = [];
    amat = sparse(arow, acol, aval, 3*(size(rbead,2)+size(rmyo,2)), 3*(size(rbead,2)+size(rmyo,2)), ind-1);
    gaval(garow==0) = [];
    garow(garow==0) = [];
    gacol(gacol==0) = [];
    gamat = sparse(garow, gacol, gaval, 3*(size(rbead,2)+size(rmyo,2)), 3*(size(rbead,2)+size(rmyo,2)), ind2-1);
    amat = amat + gamat;
end