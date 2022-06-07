fileList = {
	'kymo_prom2_1_noise.tiff',
	'kymo_prom2_2_noise.tiff',
	'kymo_prom2_3_noise.tiff',
	'kymo_prom2_4_noise.tiff',
	'kymo_prom2_5_noise.tiff',
	'kymo_prom2_6_noise.tiff',
	'kymo_prom2_7_noise.tiff',
	'kymo_prom2_8_noise.tiff',
	'kymo_prom2_9_noise.tiff',
	'kymo_prom2_10_noise.tiff'
};

N = length(fileList);
t = 1:111;
corrScale = nan(N,length(t));
c0Value = nan(N,length(t));
n = 1:N;

for l = n
    imArray2 = imread(fileList{l});
	
	imArray = imArray2;
	for k=t
		if k > size(imArray,1)
			continue
		end
		[xx,yy,c0,sigm] = corrcalcCirc((0:(length(imArray(k,:))-1)).*(0.08),movmean(double(imArray(k,:)),3));
		c0Value(l,k) = sqrt(c0)/mean(double(imArray(k,:)));
		kk = find(yy <= 0.5,1);
		corrScale(l,k) = xx(kk);
	end
end

alt = logical(ones(1,length(t)));
%alt = (t-11 < 45);

% plot fluctuation amplitude versus time
y = nanmean(c0Value);
err = nanstd(c0Value);
errorbar(t(alt)-11,y(alt),err(alt),'ro')

%{
% plot correlation length versus time
y = nanmean(corrScale);
err = nanstd(corrScale);
errorbar(t(alt)-11,y(alt),err(alt),'ro')
%}