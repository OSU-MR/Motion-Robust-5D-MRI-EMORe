function [Respiratory, Cardiac,PCIC, respRange] = extractMotionSignals(SG, opt)
%% OPTIONS
if ~isfield(opt,'nCoils');  opt.nCoils  = 12;           end     % Coil elements used for PCA (not used right now!)
if ~isfield(opt,'rLPF');    opt.rLPF    = [0.5];    end     % FIR low pass filter to select respiratory signal (Hz)
if ~isfield(opt,'cBPF');    opt.cBPF    = [0.5 3];      end     % FIR band pass filter to select cardiac signal (Hz)
if ~isfield(opt,'rcLPF');   opt.rcLPF   = [3];     end     % FIR low/band pass filter to select both respiratory and cardiac signal (Hz)
if ~isfield(opt,'nCHarm');  opt.nCHarm  = 2;            end     % # of cardiac harmonics to use (not implemented at the moment!)
if ~isfield(opt,'fOrd');    opt.fOrd    = 1001;         end     % FIR filter order (length+1), introudces transients so lower=better assuming suffiecient cuttoffs
if ~isfield(opt,'fs');      opt.fs      = 20;           end     % Self-gating sampling frequency in Hz

%% PREPROCESSING
A=zeros(size(SG));
for j = 1:size(SG,2) 
    for k = 1 : size(SG,3)        
        tmp = squeeze(SG(:,j,k))';
        % Reconstruct 1D projections of magnitude 
        A(:,j,k) = double(abs(ifftshift(ifft(fftshift(tmp),[]))));
%         A(:,j,k) = double((ifftshift(ifft(fftshift(tmp),[]))));
    end
end

% Concatenate coil elements to readout dimension
A = reshape(A,size(A,1)*size(A,2),size(A,3));

% Tranpose A for PCA -> [Time]x[Readout]
A = A'; 


% Display matrix A
figure; 
imagesc(A);
title('Casorati Matrix');
xlabel('Readout:Coils');
ylabel('Time');
savefig('Casorati.fig');
saveas(gcf,'Casorati.jpg');

% Subtract the mean from each column of A
A = A - repmat(mean(A,1),[size(A,1),1]);

% FIR filtering over temporal dimension to select respiratory/cardiac signals
rFilt = fir1(opt.fOrd-1,opt.rLPF/(opt.fs/2),hanning(opt.fOrd));     % Respiratory
cFilt = fir1(opt.fOrd-1,opt.cBPF/(opt.fs/2),hanning(opt.fOrd));     % Cardiac
% rcFilt = fir1(fOrd-1,rcLPF/(fs/2),hanning(fOrd));   % Both
rA=zeros(size(A));
cA=zeros(size(A));
for j =  1 : size(A,2)
    % Select respiratory signal
    rA(:,j) = conv(A(:,j),rFilt,'same');
    % Select cardiac signal
    cA(:,j) = conv(A(:,j),cFilt,'same');

end



%% EXTRACT RESPIRATORY AND CARDIAC MOTION USING PCA
% Principle component analysis 
% Respiratory signal
Respiratory = pcaSG(rA,1); 
respRange = 10;

nComp = 5;

% Cardiac Signal
% Take first 5 principle components
PC_C = pcaSG(cA, nComp);
% Follow up with independent component analysis
IC_C = fastICA(PC_C',nComp,'negentropy', 0); IC_C = IC_C';



% Cardiac component selection
Card_PCIC = cat(2,PC_C,IC_C);

component_rank = select_cardiac_component_ISMRM(Card_PCIC);
selected_component = component_rank(1);


% Test respiratory
PC_R = pcaSG(rA, nComp);
IC_R = fastICA(PC_R', nComp, 'negentropy',0); IC_R = IC_R';



figure;
for i = 1 : nComp
    if selected_component <= nComp && selected_component == i
        plot_color = 'r';
    else
        plot_color = 'b';
    end
    subplot(nComp,1,i);
    plot(PC_C(:,i)); title('Cardiac PCA', 'Color', plot_color);
end
savefig('CardiacPca.fig');
saveas(gcf,'CardiacPca.jpg');
figure;
for i = 1 : nComp
    if selected_component > nComp && selected_component - nComp == i
        plot_color = 'r';
    else
        plot_color = 'b';
    end
    subplot(nComp,1,i);
    plot(IC_C(:,i)); title('Cardiac ICA', 'Color', plot_color);
end
savefig('CardiacIca.fig');
saveas(gcf,'CardiacIca.jpg');


%% Extra Code to find best Cardiac Signal
PCIC = cat(2,PC_C,-PC_C,IC_C,-IC_C);
totalpeaks=zeros(20,1);
entropy=zeros(20,1);
triggers=NaN(size(PCIC));
std_rr=zeros(20,1);
edges=linspace(0.25,2,11);
figure;
for i=1:size(PCIC,2)
  tmp = -gradient(PCIC(:,i));
  tmp(tmp < 0) = 0;
  cutOff = prctile(tmp, 95)*0.5;
  tmp(tmp < cutOff) = 0;
  [pks, loc] = findpeaks(tmp);
  triggers(loc,i)=PCIC(loc,i);
  totalpeaks(i)=numel(pks);
  time_pks=loc.*opt.TR;
  rr=time_pks(2:end)-time_pks(1:end-1);
  rr(rr<0.35)=[];
  std_rr(i)=std(rr);
  hist_rr=histcounts(rr,edges);
  pdf_rr=hist_rr/numel(rr);
  entropy(i)=-sum(pdf_rr(pdf_rr>0).*log2(pdf_rr(pdf_rr>0)));
  subplot(4,5,i);
  histogram(rr,edges);
  title("E="+round(entropy(i),2)+" S="+round(std_rr(i),3));
end
savefig('hist_rr.fig');
saveas(gcf,'hist_rr.jpg');
[min_s,ind_s]=min(std_rr);
[min_e,ind_e]=min(entropy);

if(opt.liver)
Cardiac = Card_PCIC(:, component_rank(1)); 
else
Cardiac = PCIC(:,ind_e); 
end


figure;
for i = 1 : 5
    if (ind_s == i) && (ind_e==i)
        plot_color = 'r';
    elseif  (ind_s == i) || (ind_e==i)
        plot_color = 'g';
    else
        plot_color = 'b';
    end
    subplot(5,1,i);
    plot(1:length(PC_C(:,i)),PC_C(:,i),1:length(triggers(:,i)),triggers(:,i),'.r'); 
    title("PCA Entropy="+round(entropy(i),2)+" Std="+round(std_rr(i),2), 'Color', plot_color);
end
savefig('PCIC1.fig');
saveas(gcf,'PCIC1.jpg');
figure;
for i = 6 : 10
   if (ind_s == i) && (ind_e==i)
        plot_color = 'r';
    elseif  (ind_s == i) || (ind_e==i)
        plot_color = 'g';
    else
        plot_color = 'b';
    end
    subplot(5,1,i-5);
    plot(1:length(PCIC(:,i)),PCIC(:,i),1:length(triggers(:,i)),triggers(:,i),'.r'); 
    title("Inverted PCA Entropy="+round(entropy(i),2)+" Std="+round(std_rr(i),2), 'Color', plot_color);
end
savefig('PCIC2.fig');
saveas(gcf,'PCIC2.jpg');

figure;
for i = 11 : 15
   if (ind_s == i) && (ind_e==i)
        plot_color = 'r';
    elseif  (ind_s == i) || (ind_e==i)
        plot_color = 'g';
    else
        plot_color = 'b';
    end
    subplot(5,1,i-10);
    plot(1:length(PCIC(:,i)),PCIC(:,i),1:length(triggers(:,i)),triggers(:,i),'.r'); 
    title("ICA Entropy="+round(entropy(i),2)+" Std="+round(std_rr(i),2), 'Color', plot_color);
end
savefig('PCIC3.fig');
saveas(gcf,'PCIC3.jpg');

figure;
for i = 16 : 20
   if (ind_s == i) && (ind_e==i)
        plot_color = 'r';
    elseif  (ind_s == i) || (ind_e==i)
        plot_color = 'g';
    else
        plot_color = 'b';
    end
    subplot(5,1,i-15);
    plot(1:length(PCIC(:,i)),PCIC(:,i),1:length(triggers(:,i)),triggers(:,i),'.r'); 
    title("Inverted ICA Entropy="+round(entropy(i),2)+" Std="+round(std_rr(i),2), 'Color', plot_color);
end
savefig('PCIC4.fig');
saveas(gcf,'PCIC4.jpg');


%% 

%==================================================================================================
% Respiratory component selection
component_rank = 1;
selected_component = component_rank(1);

figure;
for i = 1 : nComp
    if selected_component <= nComp && selected_component == i
        plot_color = 'r';
    else
        plot_color = 'b';
    end
    subplot(nComp,1,i);
    plot(PC_R(:,i)); title('Respiratory PCA', 'Color', plot_color);
end
savefig('RespiratoryPca.fig');
saveas(gcf,'RespiratoryPca.jpg');
figure;
for i = 1 : nComp
    if selected_component > nComp && selected_component - nComp == i
        plot_color = 'r';
    else
        plot_color = 'b';
    end
    subplot(nComp,1,i);
    plot(IC_R(:,i)); title('Respiratory ICA', 'Color', plot_color);
end    
saveas(gcf,'RespiratoryIca.jpg');

return;


