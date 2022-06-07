%{
% for experimental kymos
fileList = {
	'KymoGraphD1N0000C1.tif',
	'KymoGraphD1N0003C1.tif',
	'KymoGraphD1N0006C1.tif',
	'KymoGraphD1N0006C2.tif',
	'KymoGraphD2N0003C1.tif',
	'KymoGraphD2N0003C2.tif',
	'KymoGraphD2N0004C1.tif',
	'KymoGraphD2N0006C1.tif',
	'KymoGraphD2N0006C2.tif',
	'KymoGraphD2N0006C3.tif',
};

dr = '.\';
MinProm = 50;
shift = 1;
smoothWin = 1;
%}



%{
% for simulated kymos.
fileList = {
    'kymo_prom2_6_noise.tiff',
	'kymo_prom2_7_noise.tiff',
	'kymo_prom2_8_noise.tiff',
	'kymo_prom2_9_noise.tiff',
	'kymo_prom2_10_noise.tiff',
	'kymo_prom2_26_noise.tiff',
	'kymo_prom2_27_noise.tiff',
	'kymo_prom2_28_noise.tiff',
	'kymo_prom2_29_noise.tiff',
	'kymo_prom2_30_noise.tiff'
};

dr = '.\';
MinProm = 0;
shift = 11;
smoothWin = 1;
%}


N = length(fileList);
T = 157;
cnt = nan(10,T);
for l=1:N
	[imArray,I] = imread([dr fileList{l}]);
	for k=1:min(T,size(imArray,1))
		cnt(l,k) = (size(imArray,2)-1)*0.08/sum(islocalmax(movmean(imArray(k,:),smoothWin),'MinProminence',MinProm,'MinSeparation',1));
	end
end

errorbar((1:size(cnt,2))-shift,nanmean(cnt),nanstd(cnt),'bo');
figure;
plot(cnt','bo-')