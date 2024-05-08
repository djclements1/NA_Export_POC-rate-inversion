function [XN PN C] = RTSsmoother_POC(XM,XP,PM,PP,A0);

 % Set first and last day of inversion. 
 st = min(find(~isnan(nanmean(XP,2))));
 nd = max(find(~isnan(nanmean(XP,2))));

 % Initialize the smoother matrices 
 XN = nan(size(XP));   
 PN = nan(size(PP));

 % Set the final day of the smoother to the last filter estimate
 PN(:,:,nd) = PP(:,:,nd);
 XN(nd,:) = XP(nd,:);

% This is the RTS Smoother using the same Jacobian A0 as the filter!
for t = nd-1:-1:st
    C = PP(:,:,t) * A0(:,:,t)'*PM(:,:,t+1)^-1;    % This is the smoother gain
    XN(t,:) = XP(t,:)' - C*(XM(t+1,:)' - XN(t+1,:)');   
    PN(:,:,t) = PP(:,:,t) - C*(PM(:,:,t+1) - PN(:,:,t+1))*C';
end