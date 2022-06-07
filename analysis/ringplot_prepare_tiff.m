function [] = ringplot_prepare_tiff(itask,iend)
	%{
	make a tiff file of the run with label itask
	for the runs 1:iend.	
	%}
	suffix = 'min.mat';
	indx = 1:iend;

	for n = itask;
		prefix = ['prom_' num2str(n) '_'];
		for i = indx
			iload = i;
			filename = strcat(prefix,num2str(iload),suffix);
			if exist(filename,'file')
				load(filename,'rmyo','fc','ipt','ifor','rbead','bancf','bancm','xmat');
			else
				continue;
			end
			temp1 = fluo_emu_myo(50, rmyo, bancm);
			temp1 = rgb2gray(temp1);
			if i == 1
				imwrite(temp1, ['prom2_' num2str(n) '.TIFF'])
			else
				imwrite(temp1, ['prom2_' num2str(n) '.TIFF'], 'writemode', 'append')
			end
		end
	end
end