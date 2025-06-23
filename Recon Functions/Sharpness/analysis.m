clear;
clc;

name=uipickfiles;%('filterspec', 'C:\Users\Rizwan\Research\SHARPNESS\Sharpness_Perfusion\Sharpness');
fid = fopen(char(name{1})); 
load(char(name{1}));
S=size(Sav,2);
fclose('all');

D=sparse(zeros(1e5,S));      % Concatination of all pdfs
M=zeros(size(name,2),S);     % Individual pdf means
ind=zeros(1,S);
T=zeros(S*(S-1)/2,S-1);
for i=1:size(name,2)
    fid = fopen(char(name{i})); 
    load(char(name{i}));
    for j=1:S
        D(ind(j)+1:ind(j)+numel(Sav(j).Sh),j)=Sav(j).Sh;
        M(i,j)=mean(Sav(j).Sh);
        ind(j)=ind(j)+numel(Sav(j).Sh);
    end
    T((i-1)*S+1:i*S,:)=Sav(1).t;
end
fclose('all');

%% Pair-wise t-test
t=zeros(S-1, S-1);
for k=1:S
    if k>1
        for g=1:1:k-1
            t(g,k-1)=ttest2(full(D(D(:,g)~=0)), full(D(D(:,k)~=0)),0.05,'both','unequal');
        end
    end
end
t=[t; zeros(1,k-1)]+[zeros(1,k-1); t'];


%% =========
load cmap;
c=round(linspace(1,size(cmap,1),S));
figure;
for i=1:S
    subplot(S,1,i); hist(D(D(:,i)~=0),16);
    h = findobj(gca,'Type','patch');
    set(h,'FaceColor',cmap(c(i),:), 'EdgeColor','w');
    title(['Method ', num2str(i)]);
    legend(['Mean: ' num2str(mean(D(D(:,i)~=0))) '\pm' num2str(mean(D(D(:,i)~=0))) ', H: ' num2str(t(i,:))]);
end
