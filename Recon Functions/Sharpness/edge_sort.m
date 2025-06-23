function[St] = edge_sort(I,p)

c_L=p.c_L;
P=regionprops(I,'PixelIdxList');
St(size(P,1)).SortList=randn(size(I,1),1); % Not sure if this is the best
                                             % way to initialize
I=bwlabel(I);
% figure;
for i=1:size(P,1)
    ind=P(i).PixelIdxList;
    I_tmp = I==i;
    ind_first=end_find(I_tmp); 
    if isempty(ind_first)
        loop=1;
        ind_first=P(i).PixelIdxList(1);
    elseif numel(ind_first)==2
        loop=0;
        ind_first=min(ind_first);
    else
        error('Branching detected!');
    end
    
    ind_sort=ind;
    N=numel(ind);
    ind_sort(1)=ind_first;
    ind=setdiff(ind,ind_first);
    for j=2:N
        cntr = ind_sort(j-1);
        dis=dis_pxl(cntr,ind,size(I));
        [a,b]=min(dis);
        ind_sort(j)=ind(b);
        ind(b)=[];
    end
    if loop==1
        ind_sort=[ind_sort; ind_sort(1:2*floor(c_L/2))];
    end
    St(i).SortList=ind_sort; %???
end
St=St(:);





function[dis] = dis_pxl(cntr,ind,n)

x0 = mod(cntr-1,n(1))+1;
y0 = ceil(cntr/n(1));

x = mod(ind-1,n(1))+1;
y = ceil(ind/n(1));

dis=sqrt((x-x0).^2 + (y-y0).^2);