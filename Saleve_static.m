clearvars; clearvars –global; close all;
addpath(genpath('topotoolbox-master/'),'-end');
global  parSPM;

modelInput_static() % load model parametes

% Model geometry
x=[0:parSPM.dx:parSPM.L];y=[0:parSPM.dx:parSPM.L];[X,Y]=meshgrid(x,y);siz=size(X);n=numel(X);Z=rand(siz).*1;nn=numel(X);dem=GRIDobj(X,Y,Z);iborder=find( X==min(min(X)) | X==max(max(X)) | Y==min(min(Y)) | Y==max(max(Y)) );
[ni]=directneighbours(dem);

% Maps of model parameters
parSPM.K=parSPM.K.*ones(size(X)); % Erodability map
% ind=find(X<(x(end)-x(1))/2);parSPM.K(ind)=parSPM.K(ind)./5; % If we want to test the heterogeneous K case
parSPM.P=parSPM.P.*ones(size(X)); % Precipitation map
% ind=find(X<(x(end)-x(1))/2); parSPM.P(ind)=parSPM.P(ind)./5; % If we want to test the heterogeneous P case
parSPM.U=parSPM.U.*ones(size(X)); % Uplift map
% ind=find(X<(x(end)-x(1))/2); parSPM.U(ind)=parSPM.U(ind)./5; % If we want to test the heterogeneous P case

for it=1:parSPM.Niter
    it
    %% Topotoolbox
    % Fill sinks and update water height (or discharge at outlets impose equal 0)
    dem = fillsinks(dem);
    % Compute slope
    slope = gradient8(dem);
    % Compute flow direction
    FD = FLOWobj(dem,'mex',true); % single-flow
    % Compute discharge
    A  = flowacc(FD,parSPM.P);discharge=A.Z.*parSPM.dx.^2;
    % Compute upstream distance and wave propagation time
    [updist,uptime] = flowtimedistance(FD.ix,FD.ixc,FD.size,FD.cellsize,parSPM.K,discharge,parSPM.m1,parSPM.m2,parSPM.m3,parSPM.Qc1,parSPM.Qc2);
    % Save elevation
    demold=dem; 
    %% Salève
    % River Erosion
    ixtemp  = double(FD.ix); % Donor nodes
    ixctemp = double(FD.ixc);% Receiver nodes
    start = numel(ixtemp);
    for r = start:-1:1
        dem.Z(ixtemp(r)) = dem.Z(ixctemp(r))+parSPM.U(ixctemp(r)).*(uptime(ixtemp(r))-uptime(ixctemp(r)));
    end
    % Ensure base level
    dem.Z(iborder)=0;
    % Compute crest disequilibrium
    DB = drainagebasins(FD);[crest_n,crest_dtdiff]=crestDisequilibrium(ni,dem,DB);
    % Compute erosion rate (m/d)
    E=parSPM.U.*parSPM.T-dem.Z;
    %% Plot
    if mod(it,1)==0
        subplot(1,3,1);imagesc(dem);axis square;xlabel('x (m)');ylabel('y (m)'); 
        subplot(1,3,2);loglog(discharge,slope.Z,'.k');axis square;xlabel('discharge');ylabel('slope'); 
        subplot(1,3,3);semilogx(it,nanmean(nanmean(crest_dtdiff./parSPM.dx)),'.k');hold on;axis square;xlabel('N_{iter}');ylabel('\Deltaz_{crest}');
        drawnow
    end 
end