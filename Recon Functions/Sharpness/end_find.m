function[end_ind]=end_find(I)

n=-1*ones(3,3,4*4);

n(1,2,1)=1; n(2,2,1)=1;
n(:,:,2)=imrotate(n(:,:,1),90);
n(:,:,3)=imrotate(n(:,:,1),180);
n(:,:,4)=imrotate(n(:,:,1),270);

n(1,1,5)=1; n(2,2,5)=1;
n(:,:,6)=imrotate(n(:,:,5),90);
n(:,:,7)=imrotate(n(:,:,5),180);
n(:,:,8)=imrotate(n(:,:,5),270);

n(1,1,9)=1; n(1,2,9)=1; n(2,2,9)=1;
n(:,:,10)=imrotate(n(:,:,9),90);
n(:,:,11)=imrotate(n(:,:,9),180);
n(:,:,12)=imrotate(n(:,:,9),270);

n(1,3,13)=1; n(1,2,13)=1; n(2,2,13)=1;
n(:,:,14)=imrotate(n(:,:,13),90);
n(:,:,15)=imrotate(n(:,:,13),180);
n(:,:,16)=imrotate(n(:,:,13),270);


end_ind=[];
for j=1:16
    tmp=conv2(double(I),n(:,:,j),'same');
    end_ind=[end_ind; find(tmp==sum(sum(sum( max(n(:,:,j),0)))))];
%     end_len=numel(end_ind);
end