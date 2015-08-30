for m = 2:6
  hold off;
  for basisflag=1:4
    e1(m-1,basisflag,:) = colloc(m,basisflag);
    hold on;
  end
  xlabel('$x$','Interpret','Latex','FontName','Times','FontSize',16);
  ylabel('$\epsilon$','Interpret','Latex','FontName','Times','FontSize',16);
  set(gca,'FontName','Times','FontSize',16);
  xlim=get(gca,'xLim');
  ylim=get(gca,'yLim');
  text(0.75,(ylim(2)-ylim(1))*0.1+ylim(1),strcat('$N=',num2str(m),'$'),'Interpret','Latex','FontName','Times','FontSize',16);
  legend('Chebyshev','Legendre','Lag Chebyshev','Lag Legendre','Location','NorthWest')
% pause
  savefig(strcat('colloctest2_',num2str(m)),'pdf');
end

hold off;
semilogy(e1(:,1,1),e1(:,1,2),'r',...
         e1(:,2,1),e1(:,2,2),'b',...
         e1(:,3,1),e1(:,3,2),'mx',...
         e1(:,4,1),e1(:,4,2),'cx');
xlabel('$N$','Interpret','Latex','FontName','Times','FontSize',16);
ylabel('$\|\epsilon\|$','Interpret','Latex','FontName','Times','FontSize',16);
set(gca,'FontName','Times','FontSize',16);
savefig('colloctest2conv','pdf');
