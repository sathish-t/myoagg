function [temp] = fluo_emu_myo(norm_int, rmyo, bancm)
	%{
	generate one frame of a TIFF file 
	
	arguments:
		norm_int: intensity used for normalizing image
		rmyo: myosin positions
		bancm: boolean array, 1 if respective myosin is anchored, 0 otherwise.
			please note that in this application all myosins are anchored,
			so bancm must be set to all ones.
		
	returns:
		one frame of tiff file.
	
	%}
	
	pixel_size = 0.08;
	sigma = 0.6 / 2.36; % wavelength of the microscope, used as the width of the Gaussian point-spread function
	inv_s_sq = sigma^-2;
	rmyoa = rmyo(:,bancm);
	rmyou = rmyo(:,~bancm);

	theta = linspace(0,2*pi,9);
		theta(end) = [];

	left = -5;
	right = 5;
	top = 5;
	bottom = -5;

	x_vec = left:pixel_size:right;
	y_vec = bottom:pixel_size:top;
	fluo = zeros(numel(y_vec), numel(x_vec));
	[x,y] = meshgrid(x_vec,y_vec);

	temp = zeros(size(fluo,1),size(fluo,2),3);

	fluo = zeros(numel(y_vec), numel(x_vec));
	for i = 1:size(rmyoa,2)
		this_rmyo = rmyoa(:,i);
		dx = x - this_rmyo(1);
		dy = y - this_rmyo(2);
		dr_sq = dx .* dx + dy .* dy;
		r2 = this_rmyo .* this_rmyo;
		r2 = r2(1)+r2(2);
		if r2 > 1
		fluo = fluo + 2 * exp(-dr_sq .* .5 * inv_s_sq);
		else        
		fluo = fluo + 2 * exp(-dr_sq .* .5 * inv_s_sq);
		end
	end

	temp(:,:,1) = fluo /norm_int;
	temp(:,:,2) = temp(:,:,1);
	temp(:,:,3) = temp(:,:,1);
end
