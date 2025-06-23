function card = extractCardSignal(card_caso, opt, folder_opt_path)
%EXTRACTCARDSIGNAL  Extracts & saves the cardiac trace from SG data
%
%   card = extractCardSignal(card_caso, opt, folder_opt_path)
%   Inputs:
%     card_caso       – [Time × Readouts] cardiac-filtered Casorati matrix
%     opt.TR_SG       – SG repetition time (s)
%     folder_opt_path – directory to save diagnostic figures
%
%   Output:
%     card – oriented cardiac trace (vector)

    % PCA + ICA component construction
    card_pca = pcaSG(card_caso, 5);
    card_ica = fastICA(card_pca', 5, 'negentropy', 0);
    card_ica = card_ica';
    PCIC = cat(2, card_pca, -card_pca, card_ica, -card_ica);

    % Compute and save RR-interval histograms
    totalpeaks = zeros(10,1);
    entropy   = zeros(10,1);
    triggers  = NaN(size(PCIC));
    std_rr    = zeros(10,1);
    edges     = linspace(0.25,2,11);
    hist_fig = figure('Visible','off');
    for i = 1:size(PCIC,2)
        tmp = -gradient(PCIC(:,i)); tmp(tmp<0)=0;
        cutOff = prctile(tmp,95)*0.5; tmp(tmp<cutOff)=0;
        [pks, loc] = findpeaks(tmp);
        triggers(loc,i) = PCIC(loc,i);
        totalpeaks(i)    = numel(pks);
        time_pks = loc * opt.TR_SG;
        rr       = diff(time_pks); rr(rr<0.35)=[];
        std_rr(i) = std(rr);
        hist_rr   = histcounts(rr, edges);
        pdf_rr    = hist_rr / numel(rr);
        entropy(i)= -sum(pdf_rr(pdf_rr>0) .* log2(pdf_rr(pdf_rr>0)));
        subplot(4,5,i);
        histogram(rr, edges);
        title(sprintf('E=%.2f S=%.2f', entropy(i), std_rr(i)));
    end
    close(hist_fig);

    % Plot trigger traces for PCA/ICA selections
    [~, ind_s] = min(std_rr);
    [~, ind_e] = min(entropy);
    % card1_fig = figure('Visible','off');
    % for i=1:10
    %     if (ind_s==i)&&(ind_e==i), color='r'; elseif (ind_s==i)||(ind_e==i), color='g'; else color='b'; end
    %     subplot(5,2,i);
    %     plot(1:length(PCIC(:,i)), PCIC(:,i), 1:length(triggers(:,i)), triggers(:,i), '.r');
    %     title(sprintf('PCA E=%.2f S=%.2f', entropy(i), std_rr(i)), 'Color', color);
    % end
    % close(card1_fig);
    % 
    % card2_fig = figure('Visible','off');
    % for i=1:10
    %     idx = i+10;
    %     if (ind_s==idx)&&(ind_e==idx), color='r'; elseif (ind_s==idx)||(ind_e==idx), color='g'; else color='b'; end
    %     subplot(5,2,i);
    %     plot(1:length(PCIC(:,idx)), PCIC(:,idx), 1:length(triggers(:,idx)), triggers(:,idx), '.r');
    %     title(sprintf('ICA E=%.2f S=%.2f', entropy(idx), std_rr(idx)), 'Color', color);
    % end
    % close(card2_fig);

    % Select best component by minimum entropy and orient it
    card = PCIC(:,ind_e); 
    card = orientCardiac(card, opt.TR_SG);
    disp('Cardiac signal extracted...');
end
