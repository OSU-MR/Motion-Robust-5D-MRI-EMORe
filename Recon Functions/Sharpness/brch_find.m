function[brch_ind]=brch_find(I)

% Find pixels that tie three or more edges together

% Input ==============
% I: Input binary image

% Output =============
% brch_ind: Linear index for the braching pixel


n=zeros(3,3,16); % Number of possible templates for branching

n(:,:,1)=[1 0 1; 0 1 0; 1 0 0];
n(:,:,2)=imrotate(n(:,:,1),90);
n(:,:,3)=imrotate(n(:,:,1),180);
n(:,:,4)=imrotate(n(:,:,1),270);

n(:,:,5)=[0 1 0; 0 1 1; 0 1 0];
n(:,:,6)=imrotate(n(:,:,5),90);
n(:,:,7)=imrotate(n(:,:,5),180);
n(:,:,8)=imrotate(n(:,:,5),270);

n(:,:,9)=[1 0 1; 0 1 0; 0 1 0];
n(:,:,10)=imrotate(n(:,:,9),90);
n(:,:,11)=imrotate(n(:,:,9),180);
n(:,:,12)=imrotate(n(:,:,9),270);

n(:,:,13)=[1 0 0; 0 1 1; 0 1 0];
n(:,:,14)=imrotate(n(:,:,13),90);
n(:,:,15)=imrotate(n(:,:,13),180);
n(:,:,16)=imrotate(n(:,:,13),270);


brch_ind=[];
for j=1:16
    tmp=conv2(double(I),n(:,:,j),'same');
    brch_ind=[brch_ind; find(tmp>=sum(sum(sum(n(:,:,j)))))];
%     end_len=numel(end_ind);
end