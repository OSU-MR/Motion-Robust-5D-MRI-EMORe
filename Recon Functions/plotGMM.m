function [H] = plotGMM(X, mu, sigma, phi)

h = histogram(X,100); hold on;
h.FaceAlpha = 0.2;
% h.FaceColor = 'gray';

% axis([-2,6,0,935]);
% grid on
% axis('square');

% Order X
X = sort(X,'ascend');


% X1 = linspace(min(X),coLow,100);
% X2 = linspace(coLow,coHi,100);
% X3 = linspace(coHi,max(X),100);

Z = linspace(min(X),max(X),length(X));

% % Cutoffs
% C1 = mu(2) - 3*sigma(2);
% C2 = mu(2) + 3*sigma(2);

% pdf

for i = 1 : length(mu)
G_full(:,i) = pdf('Normal',Z',mu(i),sigma(i)) * (length(X)*h.BinWidth) * phi(i);
plot(Z,G_full(:,i),'LineWidth',2.5); hold on
end
% G3_full = pdf('Normal',Z',mu(3),sigma(3)) * (length(X)*h.BinWidth) * phi(3);

% G1 = pdf('Normal',X1',mu(1),sigma(1)) * (length(X)*h.BinWidth) * phi(1);
% G2 = pdf('Normal',X2',mu(2),sigma(2)) * (length(X)*h.BinWidth) * phi(2);
% G3 = pdf('Normal',X3',mu(3),sigma(3)) * (length(X)*h.BinWidth) * phi(3);

% plot
% plot(X1,G1,'b','LineWidth',2); hold on
% plot(X2,G2,'g','LineWidth',2); hold on
% plot(X3,G3,'r','LineWidth',2); hold off

% plot(X1,G1,'b','LineWidth',2.5); hold on
% plot(X3,G3,'r','LineWidth',2.5); hold on
% plot(X2,G2,'g','LineWidth',2.5); hold on
% 
% plot(Z,G1_full,'b','LineWidth',2.5); hold on
% % plot(Z,G3_full,'r','LineWidth',2.5); hold on
% plot(Z,G2_full,'r','LineWidth',2.5); hold on

% plot([coLow,coLow],[max(h.BinCounts),max(h.BinCounts)-5],'k'); hold on;
% plot([coHi,coHi],[max(h.BinCounts),max(h.BinCounts)-5],'k'); hold off;

end

