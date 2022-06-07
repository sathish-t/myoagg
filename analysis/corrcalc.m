function [sc,c,c0,cerr] = corrcalc(s,x)
	% calculate the cross-correlation.
	
	ds = mean(diff(s));
	xbar = mean(x);
	
	N = length(x);
	c = zeros(1,N);
	cerr = zeros(1,N);
	
	currentEntries = zeros(1,N);
	
	% calculate cross-correlation.
	for k = 0:floor(N/2)
		cnt = 1;
		for t = 1:N
			nextIndx = t + k;
			if nextIndx > N
				continue;
			end
			currentEntries(t) = (x(t) - xbar) * (x(nextIndx) - xbar);
			cnt = cnt + 1;
		end
		c(k + 1) = mean(currentEntries(1:cnt-1));
		cerr(k + 1) = std(currentEntries(1:cnt-1));
	end
	
	sc = (1:length(x))-1;
	sc = sc.*ds;
	
	c0 = c(1);
	
	c = c ./ c(1);
	
end