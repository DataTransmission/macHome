  function myplotState(Var,dim1,dim2,dim3,cond)
  %function myplotState(Var1,Var2,Var3,dim,cond)
  if cond == 1
     if dim3 == []
        plot(Var(:,dim1),Var(:,dim2))
     else
        plot3(Var(:,dim1),Var(:,dim2),Var(:,dim3))
     end
  end

