lambda=4;

% Exact solution for reference
syms x real;
slambda=sqrt(lambda);
ue=sinh(slambda*(x+1))/(slambda*cosh(2*slambda));

% Some plotting points
Nplot=200;
xplot=-1:2/Nplot:1;
ueplot=subs(ue,x,xplot);

for N = 2:6
  hold off;
  for basisflag=1:2
    unplot = spectral(xplot,lambda,N,basisflag);
    err=unplot-ueplot;
    emax = max(err);
    erms = sqrt(sum(err.^2))/Nplot;
    e1(N-1,basisflag,:) =  [N erms emax];
    switch (basisflag)
      case 1
        plot(xplot,err,'r');
      case 2
        plot(xplot,err,'b');
      case 3
        plot(xplot(1:4:Nplot),err(1:4:Nplot),'mx','MarkerSize',12);
      case 4
        plot(xplot(1:4:Nplot),err(1:4:Nplot),'cx','MarkerSize',12);
    end
    hold on;
  end
  xlabel('$x$','Interpret','Latex','FontName','Times','FontSize',16);
  ylabel('$\epsilon$','Interpret','Latex','FontName','Times','FontSize',16);
  set(gca,'FontName','Times','FontSize',16);
  xlim=get(gca,'xLim');
  ylim=get(gca,'yLim');
  text(0.75,(ylim(2)-ylim(1))*0.1+ylim(1),strcat('$N=',num2str(N),'$'),'Interpret','Latex','FontName','Times','FontSize',16);
  legend('Chebyshev','Legendre','Lag Chebyshev','Lag Legendre','Location','NorthWest')
% pause
  savefig(strcat('spectraltest2_',num2str(N)),'pdf');
end

hold off;
semilogy(e1(:,1,1),e1(:,1,2),'r',...
         e1(:,2,1),e1(:,2,2),'b',...
         e1(:,3,1),e1(:,3,2),'mx',...
         e1(:,4,1),e1(:,4,2),'cx');
xlabel('$N$','Interpret','Latex','FontName','Times','FontSize',16);
ylabel('$\|\epsilon\|$','Interpret','Latex','FontName','Times','FontSize',16);
set(gca,'FontName','Times','FontSize',16);
savefig('spectral2conv','pdf');
