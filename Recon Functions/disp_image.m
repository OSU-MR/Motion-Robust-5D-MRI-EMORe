function [out] = disp_image(x,lay,title_str,gamma,clip,fignum,scale)
    x = squeeze(abs(x));

    if nargin < 2 || isempty(lay)
        cols = ceil(sqrt(size(x,3)));
        rows = ceil(size(x,3) / cols);    
    else
        rows=lay(1);
        cols=lay(2);
    end
    if nargin < 3 || isempty(title_str)
        title_str = '';
    end
    
    if nargin < 4 || isempty(gamma)
        gamma = 1;
    end
    
    if nargin < 5 || isempty(clip)
        clip = 1;
    end
    if nargin < 6 || isempty(fignum)
        out = figure;
    else
        out = figure(fignum);
    end

    if nargin < 7 || isempty(scale)
        min_scale=min(abs(x(:)));
        max_scale=max(abs(x(:))); 
    else
        min_scale=scale(1);
        max_scale=scale(2);
    end

t=tiledlayout(rows,cols);

for frame=1:size(x,3)
    nexttile;
    imagesc(abs(x(:,:,frame)).^gamma,[min_scale max_scale.*clip]);title("Frame="+frame);
    colormap("gray");
end

title(t,title_str);
end