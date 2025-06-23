function[Data] = sigmoid_fit(Data,p)

n_s = p.n_s;
n_L = p.n_L;
n_f = p.n_f;
trim = p.trim;
typ  = p.typ;

% L= (numel(St(1).SortList) - numel(St(1).Angle))/2;

R=linspace(-n_L/2,n_L/2,n_s)';% R(end)=[];

% Find vertical range of the data
vrng=zeros(size(Data,2),1);
hrng=zeros(size(Data,2),1);
for i=1:size(Data,2)
    vtmp=[];
    htmp=[];
    for j=1:numel(Data(i).Int)/n_s
        [maxv, maxl] = max(Data(i).Int((j-1)*n_s+1:j*n_s));
        [minv, minl] = min(Data(i).Int((j-1)*n_s+1:j*n_s));
        vtmp(j) = maxv-minv;
        htmp(j) = abs(maxl-minl);
    end
    vrng(i)=median(vtmp);
    hrng(i)=median(htmp);
end
%size(Data),size(Data,2), numel(Data(1).Int)/n_s, whos, pause,
for i=1:size(Data,2)
    F=[];
    D=[];
    P=[];
    for j=1:numel(Data(i).Int)/n_s
        d = Data(i).Int((j-1)*n_s+1:j*n_s);
        [dt, Rt] = exc(d,R,n_s,n_L,hrng(i),n_f,trim); %???
%         dt=d; Rt=R; %???
        if numel(dt) >  max(ceil(1*n_s/n_L), 4) ...
                    && (max(dt)-min(dt)) > 1/3*vrng(i) ...
                    && (max(dt)-min(dt)) > 1/3*(max(d)-min(d)) % Ignore short segments and not-tall-enough segments
            if strcmp(typ, 'cf')
                p0  = [0.9,     min(dt),   1.50*(max(dt)-min(dt)),     Rt(round(end/2))];
                plb = [5e-2,   -inf,       0.50*(max(dt)-min(dt)),     Rt(floor(n_s/n_L/2))];
                pub = [5e0,     inf,       3.00*(max(dt)-min(dt)),     Rt(end-floor(n_s/n_L/2))];
                options=optimset('MaxFunEvals',2e4,'MaxIter',2e3,'TolFun',1e-8,'TolX',1e-8, 'Display', 'off');
                p=lsqnonlin(@costfun,p0,plb,pub,options,Rt,dt);
                fitt=costfun(p,Rt,zeros(size(dt)));
                F=[F;fitt];
                D=[D;dt];
                P=[P;p];
            elseif strcmp(typ, 'grad')
                D=[D;dt];
                tmp=sort(abs(diff(dt))); 
                P=[P;sqrt(mean(tmp.^2))]; % ???
            elseif strcmp(typ, '1090')
                D=[D;dt];
                d10=find(dt<(min(dt)+((max(dt)-min(dt))/10)));   d10=d10(end);
                d90=find(dt>(min(dt)+(9*(max(dt)-min(dt))/10))); d90=d90(1);
                p=abs(d90-d10); 
                P=[P;1.9085*n_s/n_L/p];
            end
        end
    end
    Data(i).Param = P;
    Data(i).Fit   = F;
    Data(i).IntS  = D;
    
end
% figure; plot(D); hold on; plot(F,'r');


% ------------------------------
function [e] = costfun(p,R,d)
e = p(3)./(1+10.^(p(1)*(p(4)-R))) + p(2) - d;



% ------------------------------
function [dt, Rt] = exc(d,R,n_s,n_L,hrng,n_f,option)
if mean(d(1:floor(end/2))) > mean(d(ceil(end/2)+1:end)) % Make all profiles ascending
    d=flipud(d);
end
fl=max(2*ceil(n_s/n_L/2)+1, 2*ceil(hrng/2)+1); % filter length
% w1=gausswin(n_s,n_f); w1=w1/sum(w1(:)); % std = n_s/(2*n_f)
w2=gausswin(fl,n_f);  w2=w2/sum(w2(:));
% dc1=conv([repmat(d(1),[floor(n_s/2),1]); d; repmat(d(end),[floor(n_s/2),1])],w1,'valid');
dc2=conv([repmat(mean(d(1:1)),[floor(fl/2),1]); d; repmat(mean(d(end:end)),[floor(fl/2),1])],w2,'valid');

% [~,a1]=min(dc1); [~,b1]=max(dc1);
[~,a2]=min(dc2); [~,b2]=max(dc2);

dt=[];
Rt=[];
if a2<b2-round(1*n_s/n_L)-1
    
    if strcmp(option,'one'); % All pixels
%         Option 1
        dt=d;
        Rt=R;    
    elseif strcmp(option,'two') % All samples between max - min
%         Option 2
        dt  = d(max(a2,1) : min(b2,n_s));
        Rt  = R(max(a2,1) : min(b2,n_s));
        
    elseif strcmp(option,'three') % All samples between max - min
%         Option 3
        dt  = dc2(max(a2,1) : min(b2,n_s));
        Rt  = R(max(a2,1) : min(b2,n_s));
        
    elseif strcmp(option,'four') % All samples belonging to the largest jump
%         Option 4
        D=diff([dc2(1);dc2]);
        D(D<1e-3*max(D))=0;
        
        lab=bwlabel(D,4); 
        auc=zeros(1,max(lab));
        for j=1:max(lab)
            auc(j)=sum(D(lab==j));
        end
        [a,J]=max(auc); 
        ind=lab==J|diff([lab==J;0]);
%         ind= ind | abs(diff([ind;0])) | abs(diff([0;ind])); % Adds one point on each side
        dt=d(ind);
        Rt=R(ind);
    end
%     figure; plot(R, dc2); hold on; plot(R, d,'r'); plot(Rt,dt, 'g');
end

% tmp=-floor(n_s/2):floor((n_s-1)/2);
% w2=exp(-1/2*((tmp/(n_s/24)).^2)); 
% w2=w2/sum(w2(:));