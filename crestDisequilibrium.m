function [crest_n,crest_dtdiff]=crestDisequilibrium(ni,dem,DB)


% Find crest nodes
crest_n=zeros(DB.size);crest_dtdiff=zeros(DB.size);
for i=1:prod(DB.size)  
    wi=DB.Z(i); % Index of current watershed
    ninode=ni(i,:);ninode(ninode==0)=[]; % Index of neighbours
    wni=DB.Z(ninode); % Index of neighbour watersheds
    ninode(wni==wi)=[];
    crest_n(i)=numel(ninode);
    if isempty(ninode)==1
        crest_dtdiff(i)=nan; 
    else
        crest_dtdiff(i)=max(abs(dem.Z(i)-dem.Z(ninode)));
    end
end