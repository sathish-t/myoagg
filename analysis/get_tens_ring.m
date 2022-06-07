function tens = get_tens_ring(segtens,fitth,rbead,ifor,ipt)
% calculate total tension scalar in actin filaments accross a boundary specified
% by fitth. 

    % preallocate 
    vtens = zeros(3,length(fitth)); % vector form
    tens = zeros(1,length(fitth)); % scalar form
    % get theta for beads
    [beadth, ~,~] = cart2pol(rbead(1,:),rbead(2,:),rbead(3,:));
    % go through each segment
    for i = 1:size(rbead,2)
        % skip formins
        if any(ifor==i)
            continue
        end
        % find elements in fitth that cross this segment
            % if this segment crosses theta == pi 
            if abs(beadth(i)-beadth(i-1)) > pi
                icross = find((fitth-beadth(i)) .* (fitth-beadth(i-1)) > 0);
            % if this segment does not cross theta == pi
            else
                icross = find((fitth-beadth(i)) .* (fitth-beadth(i-1)) < 0);
            end
        % find the vector force exerted by the small theta part to the big
        % theta part
        if beadth(i) > beadth(i-1)
            force = segtens(:,i);
        else
            force = -segtens(:,i);
        end
        % if this segment crosses theta == pi, "big theta" is small theta,
        % and vice versa, so reverse the direction of force
        if abs(beadth(i)-beadth(i-1)) > pi
            force = -force;
        end
        % add this force to all fitth that cross this segment
        for j = icross
            vtens(:,j) = vtens(:,j) + force;
        end
    end
    tens = sqrt(sum(vtens.*vtens));
end