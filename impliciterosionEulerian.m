function [dem,c]=impliciterosionEulerian(X,Y,dem,FD,discharge)

global parSPM


% River Erosion
ixtemp  = double(FD.ix); % Donor nodes
ixctemp = double(FD.ixc);% Receiver nodes
start = numel(ixtemp);
for r = start:-1:1
    don=ixtemp(r);
    rec=ixctemp(r);
    dl=sqrt((X(don)-X(rec)).^2+(Y(don)-Y(rec)).^2);
    c(don)=parSPM.K(don).*discharge(don).^parSPM.m.*parSPM.dt./dl;
    dem.Z(don) = (dem.Z(don)+c(don).*dem.Z(rec))./(1+c(don));
end

% 
% K=ones(size(X)).*parSPM.K;
% Zsave=Z;c=zeros(size(Z));
% % Compute migration rate
% for k=1:numel(nstack)
%     if nstack{k}>1  
%         c(stack{k}(1))=NaN; % celerity
%         for i=2:nstack{k} 
%             % Skip the first node of a stack that should always be a
%             % self-receiver (and therefore no erosion) - this prevents from
%             % having a if condition inside the node loop
%             s=stack{k}(i);
%             r=rec(s);
%             dl=sqrt((X(s)-X(r)).^2+(Y(s)-Y(r)).^2);
%             c(s)=K(s).*discharge(s).^parSPM.m.*parSPM.dt./dl;
%             Z(s)=( Z(s) + c(s).*Z(rec(s)) )./(1+c(s));
%         end
%     end
% end
% % Don't erode hillslopes
% % Z(discharge<parSPM.Qc)=Zsave(discharge<parSPM.Qc);
% % Computer erosion
% E=Z-Zsave;