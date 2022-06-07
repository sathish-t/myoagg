function bmmatpt = ptcap(rbead, rmyo, ipt, dbead_sq, bmmat)

    %{
    edit boolean array with True at (i,j) if ith actin bead is interacting jth myosin
    to including actin pointed ends.
    
    arguments:
        rbead, rmyo - coordinates of actin beads and myosin beads
        ipt - indices of beads corresponding to actin pointed ends
        dbead_sq - square of the distance at which adjacent actin beads
            must be maintained by the constraint
        bmmat - boolean array with True at (i,j) if ith actin bead is interacting
            with jth myosin bead, excluding pointed ends.
        
    returns:
        bmmatpt: the required array.

    %}
    
    % if the bead just prior to the pointed end is within range of a myosin
    % and the pointed end is also within range of a myosin, then unbind the
    % bead, and bind the pointed end to the myosin.
    
    [row,col] = find(bmmat);
    for i = 1:numel(row)
        thisrow = row(i);
        if any(ipt == i+1)          
            thiscol = col(i);
            dr = rbead(:,thisrow) - rmyo(:,thiscol);
            if sum(dr.*dr) < dbead_sq
                row(i) = thisrow + 1;
            end
        end
    end
    bmmatpt = sparse(row,col,true(numel(row),1),size(rbead,2),size(rmyo,2));
end