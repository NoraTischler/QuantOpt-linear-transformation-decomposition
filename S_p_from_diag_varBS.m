function [ S_p ] = S_p_from_diag_varBS( D_p, pos_el,pos_ancilla,nA )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n=size(D_p,1);
S_p=eye(2*(n+nA),2*(n+nA));
el=D_p(pos_el,pos_el);

if el <=1
S_p(pos_el,pos_el)=el;
S_p(n+pos_ancilla,n+pos_ancilla)=el;
S_p(pos_el,n+pos_ancilla)=sin(acos(el));
S_p(n+pos_ancilla,pos_el)=-sin(acos(el));
S_p(pos_el+n+nA,pos_el+n+nA)=el;
S_p(2*n+nA+pos_ancilla,2*n+nA+pos_ancilla)=el;
S_p(pos_el+n+nA,2*n+nA+pos_ancilla)=sin(acos(el));
S_p(2*n+nA+pos_ancilla,pos_el+n+nA)=-sin(acos(el));
        
elseif el >1
    
S_p(pos_el,pos_el)=el;
S_p(pos_ancilla+n,pos_ancilla+n)=el;
S_p(pos_el,2*n+nA+pos_ancilla)=sinh(acosh(el));
S_p(pos_ancilla+n,pos_el+n+nA)=sinh(acosh(el));
S_p(pos_el+n+nA,pos_el+n+nA)=el;
S_p(2*n+nA+pos_ancilla,2*n+nA+pos_ancilla)=el;
S_p(pos_el+n+nA,pos_ancilla+n)=sinh(acosh(el));
S_p(2*n+nA+pos_ancilla,pos_el)=sinh(acosh(el));
        
end

end

