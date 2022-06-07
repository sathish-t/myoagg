% calculate the tension in each segment of actin filament (between beads)
% tension is a vector defined as the force on this bead exerted by the
% section close to the barbed end

% change fc into a 3-by-#seg matrix
fcmat = vec2mat(fc,3)';

% preallocate
nbead = size(rbead,2);
segtens = zeros(3,nbead+1); % + 1 for convenience

% go through each filament, from the pointed end
for i = sort(ipt,'descend')
    % go thru every bead on the filament
    for j = i:-1:i-100
        % if formin is reached, break
        if ismember(j,ifor)
            break
        end
        % calculate the segment force
        segtens(:,j) = fcmat(:,j) + segtens(:,j+1);
    end
end

% delete the additional 1 column added for convenience
segtens(:,end) = [];