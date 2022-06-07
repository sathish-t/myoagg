function rdiff = get_rdiff(rbead,~,ipt)

    %{
    return normalized vectors pointing from a bead to the next along a filament.
    for beads at pointed ends, set vector to zero.
    
    arguments:
        rbead - coordinates of actin beads
        ~ - 
        ipt - indices of beads corresponding to actin pointed ends
        
    returns:
        rdiff: the required vector.

    %}

    % vector from this bead to the next bead
    rdiff = circshift(rbead,[0,-1]) - rbead;
    % normalize 
    for i = 1:size(rdiff,2)
        rdiff(:,i) = rdiff(:,i) / norm(rdiff(:,i)); 
    end
    % do not count rdiff on ipt
    rdiff(:,ipt) = zeros(3,numel(ipt));
end