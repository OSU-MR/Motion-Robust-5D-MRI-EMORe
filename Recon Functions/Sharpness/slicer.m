function[]=slicer(I,d,cap)

% Plots every d^th slice of a 3D matrix
% Inputs:
%   I: 3D image
%   d: Every dth slice will be displayed.
%   cap: Figure caption (optional)
% ------ Rizwan Ahmad (ahmad.46@osu.edu) ---------


if nargin<2
    error('Not enough input arguments.');
end

Nz=size(I,3); % Number of slices
I=I/max(max(max(I))); % Normalize Image
a=ceil(sqrt(Nz/d)); % size of figure panel
figure;
for i=1:d:size(I,3)
    H((i-1)/d+1)=subplot(a,a,(i-1)/d+1); imagesc(I(:,:,i),[0,1]);
end
axis(H, 'equal', 'off','image');


    
if nargin==2
    axes('position',[.12 0.88 .8 .05],'Box','off','Visible','off');
    if d==1
        title(['Displaying every slice']);
    elseif d==2
        title(['Displaying every ' num2str(d) '^{nd} slice']);
    elseif d==3
        title(['Displaying every ' num2str(d) '^{rd} slice']);
    else
        title(['Displaying every ' num2str(d) '^{th} slice']);
    end
    set(get(gca,'Title'),'Visible','On');
    set(get(gca,'Title'),'FontWeight','normal');
    set(get(gca,'Title'),'FontSize',11);
    
elseif nargin==3
    axes('position',[.12 0.88 .8 .05],'Box','off','Visible','off');
    if d==1
        title([cap ', displaying every slice']);
    elseif d==2
        title([cap ', displaying every ' num2str(d) '^{nd} slice']);
    elseif d==3
        title([cap ', displaying every ' num2str(d) '^{rd} slice']);
    else
        title([cap ', displaying every ' num2str(d) '^{th} slice']);
    end
    set(get(gca,'Title'),'Visible','On');
    set(get(gca,'Title'),'FontWeight','normal');
    set(get(gca,'Title'),'FontSize',11);
    
elseif nargin>3
    error('Too many input arguments.');
end

% axes(ax);
