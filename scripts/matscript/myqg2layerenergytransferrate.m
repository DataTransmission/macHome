%before running this process.m file under your current folder, 
% 1. substitue u0pos380p10.mat by your filename.mat
%    use this command in vim   :%s/u0pos380p10.mat/filename.mat
% 2. change your TT to a maximum of 320 
% ready to run
format long;
%=====Time Space Parameter=============================
TT=37;%keep TT below 320 in case of "OUT of MEMORY"
ii=512;
jj=257;
%=============parameter===================
Rd = 25e5;
basinscl = 1800e5;
scale = basinscl/(jj-1);
Inverse_Rd2_nondim = scale^2/Rd^2;
y = 0:scale:1800e5;
U = 6 -(y-1800e5*.5).^2/380e5^2;
Uy= -2/380e5^2*(y-1800e5*.5);
H1 = 3e5; 
H2 = 3e5;
SS=(scale/Rd)^2;
S1=SS/2;
S2=SS/2;
tscale = scale/1;
dt_day = 86400/tscale;
time = 1*dt_day:8*dt_day:TT*8*dt_day;

%================x y array====================================
%y=(1:257)*1;
%x=(1:512)*1;
%[y x]=meshgrid(y,x);
%%============================================================

%=============Load psi1 psi2=================================
for t=1:TT
if t<10
    outpsi1(:,:,t)=load(strcat('000',num2str(t),'qg2_restart_psi1out.dta'));
    outpsi2(:,:,t)=load(strcat('000',num2str(t),'qg2_restart_psi2out.dta'));
elseif t<100
    outpsi1(:,:,t)=load(strcat('00',num2str(t),'qg2_restart_psi1out.dta'));
    outpsi2(:,:,t)=load(strcat('00',num2str(t),'qg2_restart_psi2out.dta'));
else
    outpsi1(:,:,t)=load(strcat('0',num2str(t),'qg2_restart_psi1out.dta'));
    outpsi2(:,:,t)=load(strcat('0',num2str(t),'qg2_restart_psi2out.dta'));
end
end
%============================================================

save('u0pos380p10.mat','outpsi1','outpsi2')


%=========================u1'==================
for t=1:TT
     for i=1:ii
         for j=2:jj-1
          u1(j,i,t)=0.5*(outpsi1(j-1,i,t)-outpsi1(j+1,i,t));
          u2(j,i,t)=0.5*(outpsi2(j-1,i,t)-outpsi2(j+1,i,t));
        end
     end
     u1(1,:,t) =0.5*(outpsi1(jj-1,:,t)-outpsi1(2,:,t));
     u1(jj,:,t)=u1(1,:,t);
end
%==============================================

save('u0pos380p10.mat','u1','u2','-append')
clear u1 u2

%=========================v1'==================
for t=1:TT
    for i=2:ii-1
         for j=1:jj
          v1(j,i,t)=0.5*(outpsi1(j,i+1,t)-outpsi1(j,i-1,t));
          v2(j,i,t)=0.5*(outpsi2(j,i+1,t)-outpsi2(j,i-1,t));
         end
     end
     v1(:,1,t) =0.5*(outpsi1(:,2,t)-outpsi1(:,ii-1,t));
     v1(:,ii,t)=v1(:,1,t);
end
%=================================================

save('u0pos380p10.mat','v1','v2','-append')
clear v1 v2
 
%======================BT STRMFNC====================
 for t=1:TT
     for i=1:ii
         for j=1:jj
          outpsib(j,i,t)=0.5*outpsi1(j,i,t)*(H1/(H1+H2))+0.5*outpsi2(j,i,t)*(H2/(H1+H2));
          end
     end
 end
%====================================================

save('u0pos380p10.mat','outpsib','-append')
%clear outpsib

%======================BC STRMFNC====================
 for t=1:TT
     for i=1:ii
         for j=1:jj
          outpsiz(j,i,t)=0.5*outpsi1(j,i,t)-0.5*outpsi2(j,i,t);
%%%          outpsiz1(j,i,t)=0.5*outpsi1(j,i,t)-0.5*outpsib(j,i,t);
%%%          outpsiz2(j,i,t)=0.5*outpsi2(j,i,t)-0.5*outpsib(j,i,t);
          end
     end
 end
%====================================================

save('u0pos380p10.mat','outpsiz*','-append')
clear outpsiz*

load('u0pos380p10.mat','u1','v1','v2')
%============Barotropic(BTE)  Energy Transform Rate==================
for t=1:TT
    BTE = 0;
    for j=1:jj
        for i=1:ii
           BTE = BTE - u1(j,i,t)*v1(j,i,t)*Uy(j);
        end
    end           
    BTEt(t) = BTE;
end

%%==the BTE below is a testing
%for t=1:TT
%    BTE = 0;
%    for j=2:jj-1
%        for i=1:ii
%	ub(j,i,t) = (-outpsib(j+1,i,t) + outpsib(j-1,i,t));
%	end
%    end
%    for j=1:jj
%        for i=2:ii-1
%	vb(j,i,t) = (outpsib(j,i+1,t) - outpsib(j,i-1,t));
%       end
%    end
%    ub(jj,:,t) = ub(jj-1,:,t);
%    ub(1,:,t)  = ub(2,:,t);
%    vb(:,ii,t) = (outpsib(:,2,t) - outpsib(:,ii-1,t));
%    vb(:,1,t)  = vb(:,ii,t);
%    for j=1:jj
%        for i=1:ii
%        BTE= BTE - ub(j,i,t)*vb(j,i,t)*Uy(j);
%        end
%    end
%    BTEt2(t) = BTE;
%end
%===============================================================================

save('u0pos380p10.mat','BTEt','-append')
clear BTEt u1 v1

load('u0pos380p10.mat','outpsiz','v*')
%==============Baroclinic(BCE) Energy Transform Rate=================================
for t=1:TT
    BCE = 0;
    for i=1:ii-1
        for j =1:jj
            BCE= BCE + 0.5*(v1(j,i,t) + v2(j,i,t))*outpsiz(j,i,t)*Inverse_Rd2_nondim*U(j);
        end
    end
    BCEt(t) = BCE;
end
%===============================================================================

save('u0pos380p10.mat','BCEt','-append')
clear outpsiz v*

%==============BTE & BCE Tnansfer Rate Time Series Plot========================================
load('u0pos380p10.mat','BCEt','BTEt')
figure
hold on;
plot(time,BTEt/((ii-1)*(jj-1)),'b'); %here we devide by Area ((ii-1)*(jj-1))
plot(time,BCEt/((ii-1)*(jj-1)),'r');
plot(time,BCEt/((ii-1)*(jj-1))+BTEt/((ii-1)*(jj-1)),'.k');
legend('BTE','BCE','TotalTransformRate');
xlabel('Nondimensionlized time')
ylabel('Nondimensionalized Energy transfer rate')
title('Theoretical QG Eddy Energy Transfer Rate (nondimensionalized)') 
print(gcf,'-djpeg','-r350','ThryEnTransferRate_u0pos380p10')
hold off;
%============== Energy transfer rate (Theoretical & Model Comparison) =================================
for t=1:TT
    epot=0;
    for j=2:jj-1
       for i=2:ii-1
       epot=epot+(0.5*outpsi1(j,i,t)-0.5*outpsi2(j,i,t))^2;
       end
    end
    epot = epot + (ii+jj-3)*(0.5*outpsi1(1,1) - 0.5*outpsi2(1,1))^2;
    epot= 0.5*SS*epot/((ii-1)*(jj-1));
    epott(t)=epot;
end
%figure
%plot(time,epott);

for t=1:TT-1
   slope(t) = (epott(t+1) - epott(t))/dt_day;
end
figure;
hold on;
plot(time(1:TT-1),slope(1:TT-1));
plot(time,BCEt/((ii-1)*(jj-1))+BTEt/((ii-1)*(jj-1)),'r');
xlabel('Nondimensionlized time')
ylabel('Nondimensionalized Energy transfer rate')
legend('Model','Theory')
title('Theoretical VS Model QG Eddy Energy Transfer Rate (nondimensionalized)') 
print(gcf,'-djpeg','-r350','ThryModelEnTransferRate_u0pos380p10')
hold off;
