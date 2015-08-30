function payment = myATT_monthly_est(totalPayment,data,nMembers,other_fees)
%function payment = myATT_monthly_est(totalPayment,data,nMembers,other_fees)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    700 min family plan charge 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
shared_total = 50; % shared by nMembers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    individual charge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
talk = 9.99; % do not discount, mothly charge for individual members
if isempty(data) % an empty data should be the default
  disp(sprintf('\n\nHsiaochu data(1)=$25, Ginos data(2)=$15, Jinghao data(3)=$0, default\n' ))
  data = [25,15,0]; % please discount
  nMembers = length(data); % Members include Hsiaochu (305-299-7934), Gino (305-299-8134), and Jinghao (305-755-7691) 
end
shared = shared_total/nMembers; % please discount. (The monthly_charge of $50 in 305-2998134 account includes monthly_charge).
if isempty(other_fees)
  other_fees = [4.11,4.11,4.11];  % do not discount, other_fees ($4.11) should be the default
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    other variable
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pay_ratio = 0.8; % paid_ratio = discounted_cost/original_cost

individualPayments = talk + (data + shared)*pay_ratio + other_fees;

%disp(['Payments for all ' num2str(nMembers) ' members are each : [' num2str(individualPayments) '] adding up to' num2str(sum(payment)) '==' num2str(totalPayment) '?'])
members = {'Hsiaochu','Gino','Jinghao'};
disp( sprintf('Payments for each member = talk + (data + shared)*pay_ratio + other_fees :\n' ))
for i=1:nMembers
  disp( sprintf('    %s = $%.2f + ($%.2f + $%.2f)*%.2f + $%.2f = $%.2f \n', members{i}, talk, data(i), shared, pay_ratio, other_fees(i), individualPayments(i) ))
end
ttlpay = sum(individualPayments);
disp( sprintf('Estimated Bill ($%.2f) == Actual Bill ($%.2f) ? ', ttlpay, totalPayment))
