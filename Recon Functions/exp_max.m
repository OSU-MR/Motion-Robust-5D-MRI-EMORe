function [nmse_hist,psnr_hist,ssim_hist,w,x] = exp_max(w,wo,y,p_em)

%% Parameters initializations

addaptive=p_em.addaptive;
m_first=p_em.m_first;
emIter=p_em.emIter;
mu=p_em.mu;
lam=p_em.lam;
mIter=p_em.mIter;
iIter=p_em.iIter;
gStp =p_em.gStp;
sig=p_em.sig;
disp_fig=p_em.disp_fig;
curr_fig = length(findobj('type','figure'));
A=p_em.A;
At=p_em.At;
tv=p_em.tv;
atv=p_em.atv;
true_images=p_em.true_images;
ind_vec=p_em.ind_vec;
numVectors=size(w,3);
screens = get(groot, 'MonitorPositions');
screen2Size = screens(2, :);
%% Vectors initializations
nmse_hist = zeros(mIter(end)*emIter-0,1);
psnr_hist = zeros(size(nmse_hist));
ssim_hist = zeros(size(nmse_hist));
x=At(w.*y);  % initial image
min_scale=min(abs(x(:)));
max_scale=max(abs(x(:)));
if(addaptive==1)
ad=(1+4*exp(-5.*(0:emIter-1)));
else
ad=ones(emIter,1);
end
d = zeros(size(x,1),size(x,2),size(x,3),3);
b = zeros(size(x,1),size(x,2),size(x,3),3);

j=0;
nmse_i=20.*log10(norm(true_images(:)-x(:))./norm(true_images(:)));
psnr_i = 20 * log10(sqrt(numel(x)).*max(abs(true_images(:))) ./ norm(true_images(:)-x(:)));
ssim_i = ssim(abs(x),abs(true_images));
if(disp_fig)
figure(curr_fig+1);
set(gcf, 'Position', screen2Size);
t=tiledlayout(3,2);
nexttile([1 2]);
hold on;
        for i = 1:numVectors
            plot(w(1,ind_vec, i), 'LineWidth', 2, 'DisplayName', ['Vector ' num2str(i)]);
        end

hold off;
legend;
title("Weighting on EM-Iter="+j);
xlabel('Index');
ylabel('Value');
for frame=1:size(x,3)
    nexttile;
    imagesc(abs(x(:,:,frame)).^1,[min_scale max_scale.*1]);title("Frame="+frame);axis("off");
    colormap("gray");
end
    title(t,"EM Recon Image on CS Iter="+j+" NMSE="+gather(nmse_i)+" PSNR="+gather(psnr_i)+" SSIM="+gather(ssim_i));
end




%% Algorithm Loop
for k = 1:emIter
%% Likelihood Step
    
    if(k>m_first)
        r=sum(abs(A(x)-y),1).^2;
        likelihood=1/(sqrt(2*pi)*sig) * exp((-r/(sig^2)));% +lam.*sum(w-wo);
        likelihood=likelihood.*wo;
        likelihood=likelihood./(sum(likelihood,3)+eps);
        w=likelihood;
 
    end

%% Maximization Step

    for i = 1:mIter(k)
        for c = 1:iIter
            gradA =  At(w.*(A(x)- y));
            gradT = mu.* atv(tv(x)-d+b);
            x = x - gStp.*(gradA + gradT);
        end
                j=j+1;
        tvx=tv(x);
        d = sth(tvx + b, lam.*ad(k)/mu);
        b = b + (tvx - d); 
        nmse_i=20.*log10(norm(true_images(:)-x(:))./norm(true_images(:)));
        psnr_i = 20 * log10(sqrt(numel(x)).*max(abs(true_images(:))) ./ norm(true_images(:)-x(:)));
        ssim_i = ssim(abs(x),abs(true_images));


        nmse_hist(j)=gather(nmse_i);
        psnr_hist(j)=gather(psnr_i);
        ssim_hist(j)=gather(ssim_i);
        if(disp_fig)
            figure(curr_fig+1);
            set(gcf, 'Position', screen2Size);
            t=tiledlayout(3,2);
            nexttile([1 2]);
            hold on;
                    for m = 1:numVectors
                        plot(w(1,ind_vec, m), 'LineWidth', 2, 'DisplayName', ['Vector ' num2str(m)]);
                    end
            
            hold off;
            legend;
            title("Weighting on EM-Iter="+k);
            xlabel('Index');
            ylabel('Value');
            for frame=1:size(x,3)
                nexttile;
                imagesc(abs(x(:,:,frame)).^1,[min_scale max_scale.*1]);title("Frame="+frame);axis("off");
                colormap("gray");
            end
            title(t,"EM Recon Image on CS Iter="+j+" NMSE="+gather(nmse_i)+" PSNR="+gather(psnr_i)+" SSIM="+gather(ssim_i));
        end
    end
end
end