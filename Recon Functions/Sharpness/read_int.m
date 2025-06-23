function[Data] = read_int(St,I)
n(1)=size(I,1);
n(2)=size(I,2);

%% =============
XI = repmat((1:n(1))',[1,n(2)]);
YI = repmat(1:n(2),   [n(1),1]);
for i=1:size(St,1)
    if ~isempty(St(i).Line)
        x=St(i).Line(:,1);
        y=St(i).Line(:,2);
    %     Data(i).Int=griddata(XI,YI,I,x,y,'linear');
        F=TriScatteredInterp(XI(:),YI(:),I(:),'natural');  
        Data(i).Int=F(x,y);
    else 
        Data(i).Int=[];
    end
end


%% ================
% for i=1:size(St,1)
%     x=St(i).Line(:,1);
%     x=x - (n(1)/2+1);
%     x=x/n(1);
%     y=St(i).Line(:,2);
%     y=y - (n(2)/2+1);
%     y=y/n(2);
%     knots=[x(:),y(:)];
%     A   = @(I) opA(I,n,knots,2,6);
%     Data(i).Int= real(A(ifftshift(fft2(I))))/n(1)/n(2);
% end
