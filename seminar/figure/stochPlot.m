figure

load('ytruth_stoch')
U1 = stoch.U;
U1b = stoch.Ub;
ytruth1 = ytruth(:,1); % 1 for X_1
%subplot(2,1,1)
figure
ylim([-0.8 0.8])
xlim([-20 20])
ylabel('U','FontSize',12,'FontWeight','bold')
xlabel('X','FontSize',12,'FontWeight','bold')
hold on
nframe= 5000;
stepsize = round(length(ytruth1)/nframe);
iframe= 1:stepsize:length(ytruth1);
%for i= 1:50
%   plot(ytruth1(iframe(i)),U1(iframe(i)),'.') %,ytruth1(iframe(i)),U1b(iframe(i)),'r')
%   h1= plot(ytruth1(iframe(i)),U1(iframe(i)),'o','MarkerFaceColor','b','MarkerSize',12)
%   filename= sprintf('stoch_cloud_%3.3d.png', i);
%   print(gcf,'-dpng',filename)%['stoch_cloud' num2str(i) '.png'])
%   delete(h1)
%end 
plot(ytruth(iframe),U1(iframe),'.')
%print(gcf,'-dpng',['stoch_cloud' num2str(i+1) '.png'])
xx=[-18:0.1:18];
U1b_plot=0.074453+0.018361.*xx-0.000253.*xx.^2+0.000002.*xx.^3;
plot(ytruth(iframe), U1b(iframe)+e_AR1(iframe,1), 'm.');hold on;
for i=1:length(xx); 
   plot(xx(i), U1b_plot(i)+e_AR1(1+(i-1)*stepsize,1), 'm.');hold on;
end
plot(xx,U1b_plot,'k')
c=[-0.4:0.05:0.4];
z=zeros(size(xx));
for i=1:length(c)
   U1b2_plot(:,i)=U1b_plot+c(i);
end
%for i =1:len 
%p=patch([X(i) X(i) X(i)-width X(i)-width],[Y(i) Y(i)+height Y(i)+height Y(i)], color(s,:)); 
%set(p,'FaceAlpha',0.125, 'EdgeColor', 'none'); 
%end
N=2;
%xx_rep=repmat(xx,N,1);
%U1_rep=repmat(U1b_plot,N,1);
%z_rep=repmat(z,N,1);
surface([xx;xx],[U1b_plot+0.4;U1b_plot-0.4],[z;z],'edgecol','no','edgealpha',.2,'facecolor','r'); alpha(0.5);
plot(xx,U1b2_plot,'r','alpha',0.5)
print(gcf,'-dpng',['stoch_cloud' num2str(i+2) '.png'])
text(1,1,'U = 0.074453 + 0.018361X - 0.000253X^2 + 0.000002X^3');
M(i+1)=getframe(gcf)
%[-2.578869742317147e-06 -2.526981707436754e-04 0.018361088020818 0.074452508136905]


%load('ytruth_stoch_10')
%U2 = stoch.U;
%U2b = stoch.Ub;
%ytruth2 = ytruth(:,1);
%subplot(2,1,2)
%plot(ytruth2,U2,'.',ytruth2,U2b,'r')
%ylabel('U','FontSize',12,'FontWeight','bold')
%xlabel('X_{a}','FontSize',12,'FontWeight','bold')
%text(1,1,'U = 0.097562 + 0.119214X + 0.002350X^2 - 0.000311X^3');
%[-3.112402976660320e-04 0.002349824438459 0.119213843847053 0.097561647871047]

%set(gcf,'color','w')
