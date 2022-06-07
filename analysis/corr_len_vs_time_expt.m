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
	'KymoGraphD2N0006C3.tif'
}

N = length(fileList);
t = 1:60;
corrScale = nan(N,length(t));
c0Value = nan(N,length(t));
n = 1:N;

for l = n
	imArray = imread(fileList{l});
	for k=t
		[xx,yy,c0,sigm] = corrcalc((0:(length(imArray(k,:))-1)).*(0.08),movmean(double(imArray(k,:)),3));
		c0Value(l,k) = sqrt(c0)/mean(double(imArray(k,:)));
		kk = find(yy <= 0.5,1);
		corrScale(l,k) = xx(kk);
	end
	l
end

errorbar(t-1,mean(corrScale),std(corrScale),'bo')
%errorbar(t-1,nanmean(c0Value),nanstd(c0Value),'ro')