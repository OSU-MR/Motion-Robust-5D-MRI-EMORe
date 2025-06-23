function [F] = edge_final(St, E)

F=zeros(size(E));
for i=1:size(St,1)
    F(St(i).fnl)=1;
end
F=logical(F);