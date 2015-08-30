function [p1 p2 PRend PRpts]=myPRPert(PR,PRrange,ndPR)
      PRend=[PR-PRrange PR+PRrange]; % end points of alpha
      p1=(PRend(2)+PRend(1))/2;  % (b + a)/2
      p2=(PRend(2)-PRend(1))/2;  % (b - a)/2
      PRpts=[PRend(1):(PRend(2)-PRend(1))/(ndPR-1):PRend(2)]'; % arbitrary chosen points of alpha (include end pts)
