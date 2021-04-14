clearvars; clearvars –global; close all;
addpath(genpath('topotoolbox-master/'),'-end');addpath(genpath('GIS_Utilities/'),'-end');
global  parSPM;

modelInput_implicit()

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

% Temporal uplift vector
%Uvec=ones(size(parSPM.t)+[0 1]).*mean(mean(parSPM.U)); Uvec(100:end)=Uvec(100:end).*2;
Uvec=ones(size(parSPM.t)+[0 1]).*mean(mean(parSPM.U));

for it=1:numel(parSPM.t)
    it
    %% Topotoolbox
    % Fill sinks and update water height (or discharge at outlets impose equal 0)
    dem = fillsinks(dem);
    % Compute slope
    slope = gradient8(dem);
    % Compute flow direction
    FD = FLOWobj(dem,'mex',true,'preprocess','none'); % single-flow
    % Compute discharge
    A  = flowacc(FD,parSPM.P);discharge=A.Z.*parSPM.dx.^2;
    % Save elevation
    demold=dem; 
    % Uplift
    dem.Z=dem.Z+parSPM.U.*parSPM.dt;dem.Z(iborder)=0;
    % River erosion
    [dem,c]=impliciterosionEulerian(X,Y,dem,FD,discharge);    
    % Ensure base level
    dem.Z(iborder)=0;
    % Compute erosion rate (m/d)
    E=(demold.Z+Uvec(it).*parSPM.dt-dem.Z)./(parSPM.dt);
    %% Plot
    if mod(it,1)==0
        subplot(1,3,1);imagesc(dem);axis square;xlabel('x (m)');ylabel('y (m)'); 
        subplot(1,3,2);loglog(discharge,slope.Z,'.k');axis square;xlabel('discharge');ylabel('slope'); 
        subplot(1,3,3);plot(parSPM.t(it)./(365*1000),mean(mean(dem.Z)),'.k');hold on;axis square;xlabel('t (kyr)');ylabel('z_{mean} (m)');
        drawnow
    end  
end