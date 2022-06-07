function [pcpq, cvec, cdot] = get_pcpq_xy(rbead,rmyo,ifor,bancf,bancm,r_ring,dbead_sq,vs_abs,dbead_first,vpol,d_for,d_myo)
    
    %{
    calculate derivatives of constraint equations i.e. partial c/partial q or pcpq
    where q are the particle coordinates and c are the constraint equations
    the constraints are:
    1. x^2 + y^2 = r^2 for anchored components i.e. restricted to surface of cylinder
    2. fixed size between adjacent beads of actin
    
    arguments:
        rbead, rmyo - actin and myosin positions
        ifor - indices of formins in rbead
        bancf, bancm - boolean array marking anchored components
        r_ring - ring radius
        dbead_sq - square of distance at which adjacent actin beads must be maintained
        vs_abs - rate at which cylinder constricts
        dbead_first - vector of distances between formins and first actin beads
        vpol - actin polymerization rate
        d_for, d_myo - distance at which formins and myosins must be maintained from the membrane.
    
    returns: 
        pcpq - derivatives of constraints
        cvec - vector of constraints (constraint equations are cvec = 0)
        cdot - vector of derivatives of c w.r.t. time
    %}

    r_for_sq = (r_ring - d_for)^2;
    r_myo_sq = (r_ring - d_myo)^2;
    
    % prepare some arrays
    
    dbead_first_sq = dbead_first .* dbead_first;
    % concatenate x and y of rbead and rmyo
    xyzcat = [rbead,rmyo];
    % concatenate bancf and bancm
    banccat = [bancf, bancm];
    % number of anchored components
    nanc = sum(banccat);
    % number of actin beads (except formins)
    nact = size(rbead,2) - numel(ifor);
    % preallocate
    pcpq1 = zeros(nanc,3*(size(rbead,2)+size(rmyo,2)));
    c1 = zeros(1,nanc);
    cd1 = c1;
    c2 = zeros(1,nact);
    cd2 = c2;
    row2 = [];
    col2 = [];
    val2 = [];

    % constraint 1: anchored particles move on a cylinder (plasma membrane)
    % indices of anchored particles
    ianc = 1:size(xyzcat,2);
    ianc = ianc(banccat);
    % form a 2 by nanc array containing 2x and 2y 
    % the factor 2 is unnecessary but kept for easier debugging
    twoxyz = 2 * xyzcat(:,ianc);
    % constraints for anchored particles
    for i = 1:nanc;
        pcpq1(i,3*(ianc(i)-1)+1) = twoxyz(1,i);
        pcpq1(i,3*(ianc(i)-1)+2) = twoxyz(2,i);
        pcpq1(i,3*(ianc(i)-1)+3) = 0*twoxyz(3,i);
        if i <= sum(bancf)
            c1(i) = xyzcat(1:2,ianc(i))' * xyzcat(1:2,ianc(i)) - r_for_sq;
        else
            c1(i) = xyzcat(1:2,ianc(i))' * xyzcat(1:2,ianc(i)) - r_myo_sq;
        end
        cd1(i) = 2 * r_ring * vs_abs;
    end
    
    %% constraint 2: distances between adjacent actin beads are fixed
    % indices of actin beads, excluding formin
    iact = 1:size(rbead,2);
    iact(ifor) = [];
    % form an array containing xy_this - xy_previous
    xyzdiff = rbead(:,:) - circshift(rbead(:,:),[0,1]);
    % exclude formin and multiply by 2
    twoxyzdiff = 2 * xyzdiff(:,iact);
    temp = sum(twoxyzdiff .* twoxyzdiff);
    c2 = .25 * temp - dbead_sq;
    % constraints for fixed actin subunit length
    for i = 1:nact        
        row2 = [row2; i;i;i;i;i;i];
        col2 = [col2; 3*(iact(i)-1)+1;3*(iact(i)-1)+2;3*iact(i);3*(iact(i)-2)+1;3*(iact(i)-2)+2;3*iact(i)-3]; 
        val2 = [val2; twoxyzdiff(1,i); twoxyzdiff(2,i);twoxyzdiff(3,i); -twoxyzdiff(1,i); -twoxyzdiff(2,i);-twoxyzdiff(3,i)];
        cd2(i) = 0;
    end
    
    % constraint is different for the first non-formin bead
    % as the distance of the first segment increases with time
    % at the rate vpol due to polymerization
    
    % indices of first non-formin actin bead
    ind = ifor + 1;
    % change to boolean
    temp = false(1,size(rbead,2));
    temp(ind) = true;
    % exclude formins in the boolean vector
    temp(ifor) = [];
    % change to indices
    ind = find(temp);
    % change c2 and cd2 for the first non-formin actin beads
    for j = 1:numel(ind)
        i = ind(j);
        c2(i) = .25 * twoxyzdiff(:,i)' * twoxyzdiff(:,i) - dbead_first_sq(j);
        cd2(i) = - 2 * dbead_first(j) * vpol;
    end
    
    % concatenate to form the final matrix
    pcpq1 = sparse(pcpq1);
    pcpq2 = sparse(row2,col2,val2,nact,3*(size(rbead,2)+size(rmyo,2)));
    pcpq = [pcpq1;pcpq2];
    cvec = [c1, c2];
    cdot = [cd1, cd2];
end