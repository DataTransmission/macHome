lambda=4;

% Exact solution for reference
syms x real;
slambda=sqrt(lambda);
ue=sinh(slambda*(x+1))/(slambda*cosh(2*slambda));

% Some plotting points
Nplot=2000;
xplot=-1:2/Nplot:1;
ueplot=subs(ue,x,xplot);

err = zeros([Nplot+1 10 4]);
e1  = zeros([9  3 4]);

hold off;
for basisflag=1:4
  for N = 2:10
    if (basisflag < 3)
%    Use Modal Expansion and Weighed inner products
      unplot = spectral2(xplot,lambda,N,basisflag);
    else
%    Use Nodal Expansion and unweighed inner products
      unplot = spectral(xplot,lambda,N,basisflag);
    end

%     Compute error measures on xplot points
    err(:,N,basisflag)=unplot-ueplot;
    emax = max(abs(err(:,N,basisflag)));
    erms = sqrt(sum(err(:,N,basisflag).^2))/Nplot;
    e1(N-1,:,basisflag) = [N erms emax];

%     Plotting Solutions
    plot(xplot,unplot,'r',xplot,ueplot,'k')
    if (N<4) 
      labeltag=strcat('$',num2str(N),'$');
      text(xplot(3*Nplot/4), unplot(3*Nplot/4),labeltag,'Interpret','Latex','FontSize',16,'HorizontalAlignment','right','Color','red')
    end
    hold on;
  end
  xlabel('$x$','Interpret','Latex','FontSize',16);
  ylabel('$u$','Interpret','Latex','FontSize',16);
  display('hit return to continue');
  pause
  if (basisflag < 3)
    savefig(strcat('modal',num2str(basisflag)),'pdf');
  else
    savefig(strcat('nodal',num2str(basisflag)),'pdf');
  end

%     Draw Convegence Curves
  hold off;
  semilogy(e1(:,1),e1(:,2),'k',e1(:,1),e1(:,3),'r')
  xlabel('$N$','Interpret','Latex','FontSize',16);
  ylabel('$\|\epsilon\|$','Interpret','Latex','Interpret','Latex','FontSize',16);
  legend('2-norm','max-norm');
  if (basisflag < 3)
    savefig(strcat('modalerr',num2str(basisflag)),'pdf');
  else
    savefig(strcat('nodalerr',num2str(basisflag)),'pdf');
  end
  display('hit return to continue');
  pause
end


for N = 2:9
  hold off;
  for basisflag=1:4
    switch (basisflag)
      case 1
        plot(xplot,err(:,N,basisflag),'r');
      case 2
        plot(xplot,err(:,N,basisflag),'b');
      case 3
        plot(xplot,err(:,N,basisflag),'m-');
        plot(xplot(1:40:Nplot),err(1:40:Nplot,N,basisflag),'mo','MarkerSize',12);
      case 4
        plot(xplot(1:40:Nplot),err(1:40:Nplot,N,basisflag),'cx','MarkerSize',12);
    end
    hold on;
  end
  xlabel('$x$','Interpret','Latex','FontSize',16);
  ylabel('$\epsilon$','Interpret','Latex','FontSize',16);
  set(gca,'FontName','Times','FontSize',16);
  xlim=get(gca,'xLim');
  ylim=get(gca,'yLim');
  text(0.75,(ylim(2)-ylim(1))*0.1+ylim(1),strcat('$N=',num2str(N),'$'),'Interpret','Latex','FontName','Times','FontSize',16);
  legend('Modal Chebyshev','Modal Legendre','Nodal Chebyshev','Nodal Legendre','Location','NorthWest')
  pause
  savefig(strcat('spectraltest2_',num2str(N)),'pdf');
end

hold off;
semilogy(e1(:,1,1),e1(:,2,1),'r',...
         e1(:,1,1),e1(:,2,2),'b',...
         e1(:,1,1),e1(:,2,3),'mo-',...
         e1(:,1,1),e1(:,2,4),'cx');
xlabel('$N$','Interpret','Latex','FontName','Times','FontSize',16);
ylabel('$\|\epsilon\|$','Interpret','Latex','FontName','Times','FontSize',16);
set(gca,'FontName','Times','FontSize',16);
savefig('spectral2conv','pdf');
