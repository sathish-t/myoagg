% plot the ring and components, to scale.
forminColor = '#0070C0';
mbColor = '#333333';

f = figure;

for cnt=1:(length(ifor))
	k1 = ifor(cnt);
	if cnt < length(ifor)
		k2 = ifor(cnt+1)-1;
	else
		k2 = length(rbead(1,:));
	end
    
    x2 = rbead(1,min(k1+1,k2));
    x1 = rbead(1,k1);
    
    y2 = rbead(2,min(k1+1,k2));
    y1 = rbead(2,k1);
    
    th2 = atan2(y2,x2);
    th1 = atan2(y1,x1);
    
    colorLeft = '#888888';
    colorRight = '#888888';
    
    
    if (sin(th2-th1)>0)    
        p1 = plot(rbead(1,k1:k2),rbead(2,k1:k2),'b-','Color',sscanf(colorLeft(2:end),'%2x%2x%2x',[1 3])/255,'LineWidth',1.5);
    else
        p1 = plot(rbead(1,k1:k2),rbead(2,k1:k2),'b-','Color',sscanf(colorRight(2:end),'%2x%2x%2x',[1 3])/255,'LineWidth',1.5);
    end
	hold on;
end


% plot myo clusters;
xLim = 4;
pxWidth = 700;
ptSize = pxWidth/(2*xLim);

scatter(rmyo(1,:),rmyo(2,:),pi*(ptSize*0.08)^2,'filled', ...
       'MarkerFaceAlpha',1/8,'MarkerFaceColor',[0.8 0.2 0]);
scatter(rbead(1,ifor),rbead(2,ifor),pi*(ptSize*0.08/3)^2,'filled', ...
       'MarkerFaceAlpha',3/8,'MarkerFaceColor',sscanf(forminColor(2:end),'%2x%2x%2x',[1 3])/255);
       
% plot circle;
th = linspace(0,2*pi,100);
plot(1.85*2*cos(th),1.85*2*sin(th),'k-','LineWidth',1.5,'Color',sscanf(mbColor(2:end),'%2x%2x%2x',[1 3])/255);
       
set(gcf,'position',[10,10,pxWidth,pxWidth]);
xlim([-xLim xLim]);
ylim([-xLim xLim]);

axis equal;
%{
axis off;
f.WindowState = 'maximized';
saveas(f,fName);
close all;
%}