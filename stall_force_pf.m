function sfpf = stall_force_pf(bmmat,fone,nsat,nhead,fhead,fheadmyp,bancm)

    %{
    calculate stall forces exerted by each myosin cluster per filament it interacts with
    
    arguments:
        bmmat - boolean matrix with True at (i,j) if ith filament is bound to jth myosin
        fone - force per filament exerted by myosin if its not saturated with actin
        nsat - threshold number of filaments at which myosin force is saturated
        nhead - number of myosin heads per myosin bead
        fhead,fheadmyp - force per head of myo2 and myp2
        bancm - boolean array with True at i if ith myosin is anchored.
    
    returns: 
        sfpf - vector of stall forces per filament
        
    NOTE: fone * nsat == nhead * fhead
    %}
    
    % number of filaments on each myosin cluster
    nfx = sum(bmmat);
    nfx = full(nfx);
    
    % for Myo2 clusters
    sfpf = (nfx <= nsat) * fone ... if at most nsat filaments on this cluster, stall force is fone
        + (nfx > nsat) * nhead * fhead ./ nfx; % if more than nsat filaments on this cluster, a total of fone*nsat is shared
    sfpf = sfpf .* bancm; % Keep myo2 entries in sfpf only; change myp2 entries in sfpf to zero.
    
    % for Myp2 clusters
    temp = (nfx <= nsat) * fone ... if at most nsat filaments on this cluster, stall force is fone
        + (nfx > nsat) * nhead * fheadmyp ./ nfx; % if more than nsat filaments on this cluster, a total of fone*nsat is shared
    sfpf = sfpf + (~bancm).*temp;   % Myp2 stall force per filament is fhead.
    
    % make a column vector
    sfpf = sfpf';
    % replace NaN's by 1e-6
    sfpf(isnan(sfpf)) = 1e-6;
    
end