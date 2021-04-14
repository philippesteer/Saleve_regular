function [ni]=directneighbours(DEM)

% Compute direct neighbours list of a regaular grid
[x,y] = ind2coord(DEM,1:numel(DEM.Z));[r,c] = coord2sub(DEM,x,y);siz=DEM.size;

ni=zeros(numel(x),8);ni_i=zeros(numel(x),8);ni_j=zeros(numel(x),8);
% Central points
ind=find(r>1 & r<siz(1) & c>1 & c<siz(2));
ni_i(ind,1:8)=[r(ind)-1  r(ind)-1  r(ind)-1  r(ind)    r(ind)    r(ind)+1  r(ind)+1  r(ind)+1];
ni_j(ind,1:8)=[c(ind)-1  c(ind)    c(ind)+1  c(ind)-1  c(ind)+1  c(ind)-1  c(ind)    c(ind)+1];
% First row (without borders)
ind=find(r==1 & c>1 & c<siz(2));      zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)    r(ind)    r(ind)+1  r(ind)+1  r(ind)+1  zer       zer       zer     ]; 
ni_j(ind,1:8)=[c(ind)-1  c(ind)+1  c(ind)-1  c(ind)    c(ind)+1  zer       zer       zer     ];
% Last row (without borders)
ind=find(r==siz(1) & c>1 & c<siz(2)); zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)-1  r(ind)-1  r(ind)-1  r(ind)    r(ind)    zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)-1  c(ind)    c(ind)+1  c(ind)-1  c(ind)+1  zer       zer       zer     ];
% First column (without borders)
ind=find(r>1 & r<siz(1) & c==1);      zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)-1  r(ind)-1  r(ind)    r(ind)+1  r(ind)+1  zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)    c(ind)+1  c(ind)+1  c(ind)    c(ind)+1  zer       zer       zer     ];
% Last column (without borders)
ind=find(r>1 & r<siz(1) & c==siz(2)); zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)-1  r(ind)-1  r(ind)    r(ind)+1  r(ind)+1  zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)-1  c(ind)    c(ind)-1  c(ind)-1  c(ind)    zer       zer       zer     ];
% Border points
ind=find(r==1 & c==1);           zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)    r(ind)+1  r(ind)+1  zer       zer       zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)+1  c(ind)    c(ind)+1  zer       zer       zer       zer       zer     ];
ind=find(r==1 & c==siz(2));      zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)    r(ind)+1  r(ind)+1  zer       zer       zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)-1  c(ind)-1  c(ind)    zer       zer       zer       zer       zer     ];
ind=find(r==siz(1) & c==1);      zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)-1  r(ind)-1  r(ind)    zer       zer       zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)    c(ind)+1  c(ind)+1  zer       zer       zer       zer       zer     ];
ind=find(r==siz(1) & c==siz(2)); zer=zeros(numel(ind),1);
ni_i(ind,1:8)=[r(ind)-1  r(ind)-1  r(ind)    zer       zer       zer       zer       zer     ];
ni_j(ind,1:8)=[c(ind)-1  c(ind)    c(ind)-1  zer       zer       zer       zer       zer     ];
% Convert to indices
ind=find(ni_i>0);ni(ind)=sub2ind(siz,ni_i(ind),ni_j(ind));
