function bmmat = bm_nbr (rbead, rmyo, ipt, rcap)

    %{
    return boolean array with True at (i,j) if ith actin bead is interacting jth myosin,
    but exclude pointed ends.
    
    arguments:
        rbead, rmyo - coordinates of actin beads and myosin beads
        ipt - indices of beads corresponding to actin pointed ends
        rcap - myosin capture radius. 
        
    returns:
        bmmat: the required array.

    %}

    % find out possible neighbors whose x and y positions are within rcap
    [xb, xm] = meshgrid(rbead(1,:), rmyo(1,:));
    [yb, ym] = meshgrid(rbead(2,:),rmyo(2,:));    
    [zb, zm] = meshgrid(rbead(3,:),rmyo(3,:));
    xd = xb - xm;
    yd = yb - ym;
    zd = zb - zm;
    rdsq = xd .* xd + yd .* yd + zd .* zd;
    temp = rdsq < rcap * rcap;
    bmmat = temp';
    
   
    % if more than one beads on the same filament are interacting with the same
    % myosin cluster, only the one that is closer to the pointed end but is
    % not the pointed itself counts (for simplicity)
    bmmat(ipt,:) = false;
    temp = ~circshift(bmmat,[-1,0]);
    bmmat = and(bmmat, temp);
    
end