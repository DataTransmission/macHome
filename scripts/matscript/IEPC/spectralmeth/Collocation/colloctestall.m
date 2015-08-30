for basisflag=1:4
  hold off;
  for m = 2:6
    e1(m-1,:) = colloc(m,basisflag)
    hold on;
  end
  xlabel('$x$','Interpret','Latex','FontSize',16);
  ylabel('$u$','Interpret','Latex','FontSize',16);
  pause
  savefig(strcat('colloc',num2str(basisflag)),'pdf');
  hold off;
  semilogy(e1(:,1),e1(:,2),'k',e1(:,1),e1(:,3),'r')
  xlabel('$N$','Interpret','Latex','FontSize',16);
  ylabel('$\|\epsilon\|$','Interpret','Latex','Interpret','Latex','FontSize',16);
  legend('2-norm','max-norm');
  savefig(strcat('collocerr',num2str(basisflag)),'pdf');
  pause
end
