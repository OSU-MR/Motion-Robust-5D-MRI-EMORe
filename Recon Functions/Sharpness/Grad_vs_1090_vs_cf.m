clear;
clc;
close all;

sig=[0,1,2,4,8,16,32,64,128,256]*1e-3;

% 1090 (mean,median,std), grad (mean, median,std), cf (mean,median,std)

X=[0.3784	0.3784	0	    1.5541	1.5541	0	    1.5386	1.5386	0
   0.3783	0.3783	0.0006	1.5541	1.5541	0	    1.5389	1.5391	0.0052
   0.3784	0.3785	0.0012	1.5518	1.5541	0.0179	1.5389	1.5385	0.0096
   0.3783	0.3783	0.0021	1.5286	1.5541	0.0804	1.5369	1.5389	0.022
   0.3781	0.3778	0.0048	1.4783	1.4128	0.1324	1.5404	1.5394	0.0399
   0.3786	0.3784	0.0088	1.4058	1.3728	0.1972	1.5317	1.5292	0.0788
   0.3827	0.3816	0.0169	1.2524	1.1954	0.2714	1.5438	1.5382	0.1585
   0.3911	0.3908	0.0342	1.1151	1.073	0.3715	1.5685	1.5374	0.3157
   0.4214	0.4104	0.0664	0.9299	0.8179	0.3909	1.5823	1.5163	0.5734
   0.5273	0.5167	0.1253	0.9016	0.7585	0.4695	1.7043	1.5441	0.9752];

Y=X; 

Y(:,1)=X(:,1)*X(1,8)/X(1,1);
Y(:,2)=X(:,2)*X(1,8)/X(1,2);
Y(:,3)=X(:,3)*X(1,8)/X(1,2);
 
Y(:,4)=X(:,4)*X(1,8)/X(1,4);
Y(:,5)=X(:,5)*X(1,8)/X(1,5);
Y(:,6)=X(:,6)*X(1,8)/X(1,5);

Y(:,3)=Y(:,3)/4;
Y(:,6)=Y(:,6)/4; 
Y(:,9)=Y(:,9)/6;


%% mean pm std
% errorbar(sig,Y(:,1),Y(:,3),'b'); hold on; 
% errorbar(sig,Y(:,4),Y(:,6),'r'); 
% errorbar(sig,Y(:,7),Y(:,9),'g');
% errorbarlogx(10); axis([0.5e-3 256e-3 0.4 2.2]);


%% median pm std
% bar([Y(:,2),Y(:,5),Y(:,8)]);  
errorbar(Y(:,2),Y(:,3),'sb','MarkerFaceColor','blue'); hold on; 
errorbar(Y(:,5),Y(:,6),'^r','MarkerFaceColor','red'); 
errorbar(Y(:,8),Y(:,9),'vg','MarkerFaceColor','green'); axis([0 11 0.4 2.4])
% errorbarlogx(10); axis([0.5e-3 512e-3 0.4 2.4]);




% sig=[0,1,2,4,8,16,32,64,128]*1e-3;
% 
% X=[0.0172	0        1.6449	0        1.5199   0
%    0.0172	0        1.622	0.0474	 1.5218   0.0051
%    0.0172	0.0001	 1.5466	0.0749	 1.5211	  0.0092
%    0.0172	0.0001	 1.4996	0.0726	 1.5200	  0.02
%    0.0173	0.0002	 1.4303	0.1210	 1.5147	  0.0372
%    0.0173	0.0004	 1.3125	0.1825	 1.5094	  0.0741
%    0.0176	0.0008	 1.1427	0.2692	 1.5262	  0.1563
%    0.0184	0.0016	 0.8643	0.3338	 1.5024	  0.2839
%    0.0217	0.0032	 0.7459	0.4274	 1.4398	  0.4139];
% 
% X(:,2)=X(:,2)*X(1,5)/X(1,1);
% X(:,1)=X(:,1)*X(1,5)/X(1,1); 
% X(:,4)=X(:,4)*X(1,5)/X(1,3);
% X(:,3)=X(:,3)*X(1,5)/X(1,3); 
% X(:,2)=X(:,2)/2;
% X(:,4)=X(:,4)/2; 
% X(:,6)=X(:,6)/2;
% 
% errorbar(sig,X(:,1),X(:,2)); hold on; 
% errorbar(sig,X(:,3),X(:,4),'r'); 
% errorbar(sig,X(:,5),X(:,6),'g');
% errorbarlogx(10); axis([0.5e-3 256e-3 0.4 2.2]);