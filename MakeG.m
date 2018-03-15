function [ G ] = MakeG( Gsize )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
G=eye(Gsize);
for n=Gsize/2+1:Gsize
    G(n,n)=-1;
end

end

