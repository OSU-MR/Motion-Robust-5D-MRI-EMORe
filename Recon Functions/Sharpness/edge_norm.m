function[St] = edge_norm(St,Ic,p,typ)
% Find normals to each pixel on the edge
% typ-B: From the edge normals
% typ-A: From the image gradient


% Check format
if nargin==2
    c_Lb2=1; % Default
    n_L=5;
    n_s=8*n_L;
    Rng=[-pi/2, pi/2];
    typ='A';
elseif nargin==3
    typ = 'A';
elseif nargin==3 || nargin==4
    n_L = p.n_L;
    n_s = p.n_s;
    c_L = p.c_L;
    Rng = p.AngR;
    if c_L<3
        error('Parameter "len" can only assume odd values > 1');
    elseif mod(n_L,2)==0
        c_Lb2=floor((c_L-1)/2);
    else
        c_Lb2=floor(c_L/2);
    end
elseif nargin>4
    error('Incorrect number of input arguments');
end

n=size(Ic);
R=linspace(-n_L/2,n_L/2,n_s)'; %R(end)=[];

if strcmp(typ,'B')
    for i=1:size(St,1)
        X=[]; %???
        Y=[]; %???
        loc=[];
    %     Y0=[];
        for j=1+c_Lb2:numel(St(i).SortList)-c_Lb2
            xx = mod(St(i).SortList(j-c_Lb2:j+c_Lb2),n(1))+1;
    %         xx(9)=xx(5); xx(8)=xx(4); xx(7)=xx(4); xx(6)=xx(3); xx(5)=xx(3); xx(4)=xx(3); xx(3)=xx(2);
            yy = ceil(St(i).SortList(j-c_Lb2:j+c_Lb2)/n(1));
    %         yy(9)=yy(5); yy(8)=yy(4); yy(7)=yy(4); yy(6)=yy(3); yy(5)=yy(3); yy(4)=yy(3); yy(3)=yy(2);
            [err,P] = tls(xx,yy,'no'); % Total least-squares
            St(i).Angle(j-c_Lb2) = P(1); %???
            for k=1:size(Rng,1)
                if (P(1)-pi/2)>=Rng(k,1) && (P(1)-pi/2)<=Rng(k,2)
                    x0 = mod(St(i).SortList(j),n(1))+1;
                    y0 = ceil(St(i).SortList(j)/n(1));
                    x = x0 + R*cos(P(1));
                    y = y0 + R*sin(P(1));
                    loc=[loc;St(i).SortList(j)];
            %             Y0=[Y0;y0];
                    X=[X;x];
                    Y=[Y;y];
                end
            end
        end
        St(i).Line = [X,Y]; % Co-ordinates of normal lines
        St(i).fnl  = loc;
    end
    
elseif strcmp(typ,'A')
    % Find gradient map
    [gy,gx]=gradient(real(ifft2(fft2(Ic).*fft2(fftshift(gauss2d(size(Ic),p.e_sig))))));
    A=-atan(gx./gy);

    for i=1:size(St,1)
        X=[]; %???
        Y=[]; %???
        loc=[];
    %     Y0=[];
        for j=1:numel(St(i).SortList)
    %         xx = mod(St(i).SortList(j-c_Lb2:j+c_Lb2),n(1))+1;
    %         xx(9)=xx(5); xx(8)=xx(4); xx(7)=xx(4); xx(6)=xx(3); xx(5)=xx(3); xx(4)=xx(3); xx(3)=xx(2);
    %         yy = ceil(St(i).SortList(j-c_Lb2:j+c_Lb2)/n(1));
    %         yy(9)=yy(5); yy(8)=yy(4); yy(7)=yy(4); yy(6)=yy(3); yy(5)=yy(3); yy(4)=yy(3); yy(3)=yy(2);
    %         [err,P] = tls(xx,yy,'no'); % Total least-squares
            P=A(St(i).SortList(j));
            St(i).Angle(j) = P; %???
            for k=1:size(Rng,1)
                if (P)>=Rng(k,1) && (P)<=Rng(k,2)
                x0 = mod(St(i).SortList(j),n(1))+1;
                y0 = ceil(St(i).SortList(j)/n(1));
                x = x0 + R*cos(P+pi/2);
                y = y0 + R*sin(P+pi/2);
                loc=[loc;St(i).SortList(j)];
    %             Y0=[Y0;y0];
                X=[X;x];
                Y=[Y;y];
                end
            end
        end
        St(i).Line = [X,Y];
        St(i).fnl  = loc;
    end
    
else
    error('Only A and B are acceptable options');
end

