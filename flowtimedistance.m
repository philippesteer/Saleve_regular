function [D,T] = flowtimedistance(ix,ixc,siz,cellsize,K1,A,m1,m2,m3,Ac1,Ac2) 

%FLOWDISTANCE flow distance in upstream and downstream direction
%
% Syntax
%
%     D = flowdistance(FD)
%     D = flowdistance(FD,GRIDobj)
%     D = flowdistance(FD,S)
%     D = flowdistance(FD,IX)
%     D = flowdistance(FD,x,y)
%     D = flowdistance(...,direction)
%
% Description
%
%     flowdistance calculates the horizontal distance along the flow
%     network in either upstream (default) or downstream direction.
%     Downstream flow distance is the maximum distance along the drainage
%     network (flow length). If only a flow direction object (FLOWobj) is
%     supplied to the function, flowdistance calculates the distance from
%     outlets and ridges in upstream and downstream direction,
%     respectively. If seed locations for the distance transform are
%     supplied (e.g. logical GRIDobj, STREAMobj, linear index or coordinate 
%     pairs) the distance is calculated from the seed locations. When using
%     a STREAMobj as input, usually upstream distance makes sense only.
%
% Input arguments
%
%     FD          flow direction (FLOWobj)
%     GRIDobj     logical grid (GRIDobj)
%     S           instance of STREAMobj
%     IX          linear index of seeds
%     x,y         x- and y-coordinate vectors
%     direction   'upstream' (default) or 'downstream' (maximum flow path
%                 length)
%
% Output arguments
%
%     D           distance grid (GRIDobj)
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     D = flowdistance(FD);
%     imageschs(DEM,D)
%
% See also: FLOWobj, FLOWobj/flowacc, GRIDobj
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 4. March, 2016


% 4/3/2016: the function now makes copies of FD.ix and FD.ixc (see 
% FLOWobj/flowacc

%% Do calculation
% Compute distance
d = ones(size(ix),'double')*cast(cellsize,'double');
diagdistance = cast(norm([cellsize cellsize]),'double');
nrrows = siz(1);
if mod(nrrows,2) == 1
    % odd number of rows
    d(mod(ix,2) == mod(ixc,2)) = diagdistance;
else
    % even number of rows 
    d(mod(ix,2) ~= mod(ixc,2) & ~(imabsdiff(ix,ixc)==1)) = diagdistance; % ~(ix-ixc == 1 | ixc-ix == 1))
end
DIST=d;
% Compute cumulated upstream distance and time
cl   = class(DIST);
D    = zeros(siz,cl);
T    = zeros(siz,cl);
start = numel(ix); K1save=K1(ix(start));K2save=K1(ix(start)).*Ac1.^(m1-m2);
for r = start:-1:1
    dlength=DIST(r);
    % Compute propagation time increment
    % River stream power
    if A(ix(r))>Ac1
        % River stream power
        dtime=dlength./(K1(ix(r)).*A(ix(r)).^m1);
        K1save=K1(ix(r)); % save last value of K when going upstream
    elseif A(ix(r))>Ac2
        % Colluvial stream power
        K2=K1save.*Ac1.^(m1-m2); % Compute K2 as a function of last value of K1 and Ac
        dtime=dlength./(K2.*A(ix(r)).^m2);
        K2save=K2; % save last value of K when going upstream
    else
        % Colluvial stream power
        K3=K2save.*Ac2.^(m2-m3); % Compute K2 as a function of last value of K1 and Ac
        dtime=dlength./(K3.*A(ix(r)).^m3);        
    end
    D(ix(r)) = D(ixc(r))+dlength;
    T(ix(r)) = T(ixc(r))+dtime;
end