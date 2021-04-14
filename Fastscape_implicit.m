clearvars; clearvars –global; close all;
addpath(genpath('topotoolbox-master/'),'-end');

global  parSPM;

modelInput()

% Model geometry
x=[0:parSPM.dx:parSPM.L];y=[0:parSPM.dx:parSPM.L];[X,Y]=meshgrid(x,y);siz=size(X);n=numel(X);Z=rand(siz).*0.1;nn=numel(X);
DEM=GRIDobj(X,Y,Z);
[ni]=directneighbours(DEM);
% for k=1:nn
%     [i,j]=ind2sub(siz,k);
%     temp_sub=[i+1,j+1;i+1,j;i+1,j-1;i,j+1;i,j-1;i-1,j+1;i-1,j;i-1,j-1];[iremove,jremove]=find(temp_sub(:,1)<1 | temp_sub(:,1)>siz(1) | temp_sub(:,2)<1 | temp_sub(:,2)>siz(1));
%     temp_sub(unique(iremove),:)=[];
%     ni{k}=sub2ind(siz,temp_sub(:,1),temp_sub(:,2));
% end

surface=parSPM.dx.^2.*ones(size(X));
mask=ones(size(X));

% Initialize variable
Usave=zeros(size(parSPM.t));Ersave=zeros(size(parSPM.t));Ehsave=zeros(size(parSPM.t));Zmean=zeros(size(parSPM.t));Zmin=zeros(size(parSPM.t));Zmax=zeros(size(parSPM.t));CFLmean=zeros(size(parSPM.t));CFLmin=zeros(size(parSPM.t));CFLmax=zeros(size(parSPM.t));
itsave=0;
for it=1:parSPM.nt
    % Compute greatest slope
    [slope,ASP] = gradient8(Z,parSPM.dx);
    % Update receiver list
    [rec,donor,ndon,ioutlet,ioutlet_logical]=receiverdonorssinks(siz,ASP);
    % Build Fastscape stack
    [stack,rstack,nstack,label]=buildstack(siz,Z,ioutlet,ioutlet_logical,ndon,donor);
    % Compute drainage and discharge
    P=ones(siz).*parSPM.P;
    [area,discharge]=drainageArea(siz,rec,rstack,nstack);
    % Uplift
    U=ones(siz).*parSPM.U; % m/yr;
    Z(2:end-1,2:end-1)=Z(2:end-1,2:end-1)+U(2:end-1,2:end-1).*parSPM.dt;
    % River Erosion
    [Z,Er,CFL]=impliciterosionEulerian(X,Y,Z,discharge,stack,nstack,rec);  
    % Hillslope avalanche
    [Z,Eh]=hillslopeavalanche(Z);  
    %[Z,Eh,landslide]=erosionLandslide(nn,X,Y,Z,ni,surface,mask);
    
    % CFL
    if mod(it,100)==0
        itsave=itsave+1;
        imagesc(Z); drawnow
        Zsave{itsave}=Z;
    end
    Usave(it)   = mean(mean(U(2:end-1,2:end-1))).*parSPM.dt;
    Ersave(it)  = mean(mean(Er(2:end-1,2:end-1)));
    %Ehsave(it)  = mean(mean(Eh(2:end-1,2:end-1)));   
    Zmean(it)   = mean(mean(Z));
    Zmax(it)    = max(max(Z));
    Zmin(it)    = min(min(Z));
    CFLmean(it) = mean(mean(CFL));
    CFLmax(it)  = max(max(CFL));
    CFLmin(it)  = min(min(CFL));    
end

fig = figure;Zmax= max(max(Zsave{itsave}));
for idx = 1:itsave
    imagesc(Zsave{idx});colormap('parula');caxis([0 Zmax]);axis off
    drawnow
    frame = getframe(fig);
    im{idx} = frame2im(frame);
end
close;
filename = 'testAnimated.gif'; % Specify the output file name
for idx = 1:itsave
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end

% 
k=365.*1000./parSPM.dt;
subplot(2,3,[1 2]);plot(parSPM.t./(365*1000),Usave.*k,'-k',parSPM.t./(365*1000),-Ersave.*k,'-b',parSPM.t./(365*1000),-Ehsave.*k,'-b',parSPM.t./(365*1000),-(Ersave+Ehsave).*k,'-b');xlabel('t (kan)');ylabel('U ou E (mm/yr)');set(gca,'fontsize', 8)
nom=['Erosion_uplift_time'];print('-dpdf','-painters',nom);close;
% 
subplot(2,3,[1 2]);plot(parSPM.t./(365*1000),Zmean,'-k',parSPM.t./(365*1000),Zmax,'-k');xlabel('t (kan)');ylabel('Z (m)');set(gca,'fontsize', 8)
% nom=['Z_time'];print('-dpdf','-painters',nom);close;
% 
% 
% subplot(3,3,[1 2 5 6]);imagesc(Zsave{10}); axis equal tight off;colorbar;caxis([0 1000]);ylabel(colorbar, 'Z (m)');
% nom=['Zmap_10'];print('-dpdf','-painters',nom);close;
% 
% subplot(3,3,[1 2 5 6]);imagesc(Zsave{20}); axis equal tight off;colorbar;caxis([0 1000]);ylabel(colorbar, 'Z (m)');
% nom=['Zmap_20'];print('-dpdf','-painters',nom);close;
% 
% subplot(3,3,[1 2 5 6]);imagesc(Zsave{30}); axis equal tight off;colorbar;caxis([0 1000]);ylabel(colorbar, 'Z (m)');
% nom=['Zmap_30'];print('-dpdf','-painters',nom);close;
% 
% subplot(3,3,[1 2 5 6]);imagesc(Zsave{100}); axis equal tight off;colorbar;caxis([0 1000]);ylabel(colorbar, 'Z (m)');
% nom=['Zmap_100'];print('-dpdf','-painters',nom);close;
% 
% 
% 
% imagesc(Zsave{100})
% 
% 
