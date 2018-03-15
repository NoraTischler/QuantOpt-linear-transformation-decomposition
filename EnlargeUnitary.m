function [ S ] = EnlargeUnitary( U , nA )
%ENLARGEUNITARY takes a unitary matrix and provides the quasiunitary
%matrix, with a specified number of ancilla modes included
%   Inputs: 
%   U, the unitary matrix to be enlarged
%   nA, the number of ancilla modes to be inluded
%   Output:
%   S, the quasiunitary matrix

Usize=size(U,1);
Ssize=(Usize+nA)*2;
S=eye(Ssize);

S(1:Usize,1:Usize)=U;
S(Ssize/2+1:Ssize/2+Usize,Ssize/2+1:Ssize/2+Usize)=conj(U);

%test quasiunitarity of S
G=MakeG( Ssize );

if norm(S*G*S'-G)>1E-8
    error('The S from EnlargeUnitary is not quasiunitary')
end

end

