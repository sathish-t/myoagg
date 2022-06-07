%% Calculate tension
get_segtens                     % calculate the tension in each segment of actin filament (between beads)
fitth = linspace(-pi,pi,101);
fitth(end) = [];                % angular positions where tension is calculated
tens = get_tens_ring(segtens,fitth,rbead,ifor,ipt);   % calculate ring tension accross angle fitth