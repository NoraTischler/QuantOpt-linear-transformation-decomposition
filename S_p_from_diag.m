function [ S_p ] = S_p_from_diag( D_p , pos_el , pos_ancilla , nA , varargin )
%S_P_FROM_DIAG creates a quasiunitary transformation matrix corresponding
%to the single-mode modulation corresponding to ba specified element of a
%diagonal matrix
%   Inputs: 
%   D_p, a diagonal matrix
%   pos_el, the mode of interest, for which the modulation is imlpemented
%   pos_ancilla, the mode that plays the role of the ancilla
%   nA, the total number of ancilla modes
%   varargin: optional additional inputs are whether the allowed type of
%   beam splitter is 'variable' or 'balanced'
%   Output:
%   S_p, the quasiunitary matrix of the single-mode modulation

%% PRELIMINARIES
%Default choice for the type of allowed beam splitters, if not specified
if nargin==1
    BStype='variable';
else
    BStype=varargin{1};
end

n=size(D_p,1);
el=D_p(pos_el,pos_el); %the element that specifies the modulation
num_operations=1;
m=1;
if (strcmp(BStype,'balanced')&& el <1)
    num_operations=6;
end
S_p=repmat(eye(2*(n+nA),2*(n+nA)),[1,1,num_operations]);

%% CONSTRUCTION OF TRANSFORMATION MATRIX
if el <1 %loss, so a beam splitter is needed
    switch BStype
        case 'variable' %if variable beam splitters are allowed, only one matrix is required
            S_p(pos_el,pos_el)=el;
            S_p(n+pos_ancilla,n+pos_ancilla)=el;
            S_p(pos_el,n+pos_ancilla)=sin(acos(el));
            S_p(n+pos_ancilla,pos_el)=-sin(acos(el));
            S_p(pos_el+n+nA,pos_el+n+nA)=el;
            S_p(2*n+nA+pos_ancilla,2*n+nA+pos_ancilla)=el;
            S_p(pos_el+n+nA,2*n+nA+pos_ancilla)=sin(acos(el));
            S_p(2*n+nA+pos_ancilla,pos_el+n+nA)=-sin(acos(el));
        case 'balanced' %if only 50-50 beam splitters are allowed, a Mach-Zehnder interferometer is required to implement the modulation
            S_p(pos_el,pos_el,m)=exp(1i*pi/2);
            S_p(pos_el+n+nA,pos_el+n+nA,m)=exp(-1i*pi/2);
            m=m+1;
            S_p(pos_el,pos_el,m)=cos(pi/4);
            S_p(n+pos_ancilla,n+pos_ancilla,m)=cos(pi/4);
            S_p(pos_el,n+pos_ancilla,m)=sin(pi/4);
            S_p(n+pos_ancilla,pos_el,m)=-sin(pi/4);
            S_p(pos_el+n+nA,pos_el+n+nA,m)=cos(pi/4);
            S_p(2*n+nA+pos_ancilla,2*n+nA+pos_ancilla,m)=cos(pi/4);
            S_p(pos_el+n+nA,2*n+nA+pos_ancilla,m)=sin(pi/4);
            S_p(2*n+nA+pos_ancilla,pos_el+n+nA,m)=-sin(pi/4);
            m=m+1;
            S_p(pos_el,pos_el,m)=exp(1i*(acos(el)+pi));
            S_p(pos_el+n+nA,pos_el+n+nA,m)=exp(-1i*(acos(el)+pi));
            m=m+1;
            S_p(n+pos_ancilla,n+pos_ancilla,m)=exp(1i*(-acos(el)));
            S_p(2*n+nA+pos_ancilla,2*n+nA+pos_ancilla,m)=exp(-1i*(-acos(el)));
            m=m+1;
            S_p(pos_el,pos_el,m)=cos(pi/4);
            S_p(n+pos_ancilla,n+pos_ancilla,m)=cos(pi/4);
            S_p(pos_el,n+pos_ancilla,m)=sin(pi/4);
            S_p(n+pos_ancilla,pos_el,m)=-sin(pi/4);
            S_p(pos_el+n+nA,pos_el+n+nA,m)=cos(pi/4);
            S_p(2*n+nA+pos_ancilla,2*n+nA+pos_ancilla,m)=cos(pi/4);
            S_p(pos_el+n+nA,2*n+nA+pos_ancilla,m)=sin(pi/4);
            S_p(2*n+nA+pos_ancilla,pos_el+n+nA,m)=-sin(pi/4);
            m=m+1;
            S_p(pos_el,pos_el,m)=exp(1i*(pi/2));
            S_p(pos_el+n+nA,pos_el+n+nA,m)=exp(-1i*(pi/2));
    end
    
elseif el >=1 %gain, so a parametric amplifier is needed
    
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

