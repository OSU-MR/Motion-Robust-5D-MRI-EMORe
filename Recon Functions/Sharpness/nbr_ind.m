function[N] = nbr_ind(I,ind)

% Finds indices of non-zero neighbors in the 8-neighborhood of the 
% pixel with index ind

n=size(I);

N(1) = ind - (n(1)+1);
N(2) = ind -  n(1);
N(3) = ind - (n(1)-1);
N(4) = ind - 1;
N(5) = ind + 1;
N(6) = ind + (n(1)-1);
N(7) = ind +  n(1);
N(8) = ind + (n(1)+1);
N(9) = ind;

N(N<1)=[];
N(N>n(1)*n(2))=[];


N(I(N)==0)=[];