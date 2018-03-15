function [ Transformation_matrix_stack ] = ReckDecomposition_varBS( U )
%RECKDECOMPOSITION_VARBS Provides the decomposition of a unitary matrix (U) into
%a sequence of basic 2x2 actions as outlined in the PRL
%https://doi.org/10.1103/PhysRevLett.73.58
%   The basic 2x2 actions are variable beam splitters (BS) and phase shifters
%   (PS). The output is given in terms of matrices that
%   describe the actions, in the order to be applied according to the third
%   dimension of the 3D array.

%% PRELIMINARIES AND INITIALISATION
% Check unitarity
M=size(U,1);
if norm(U*U'-eye(M))+norm(U'*U-eye(M))>1E-8
    error('input into ReckDecomposition is not unitary')
end

%initialise variables
U_current=U;
phis=NaN*ones(M,M);
thetas=NaN*ones(M,M);
m=1;
%Transformation_summary = {'transformation number','transformation', 'first mode', 'second mode'};

%% DECOMPOSITION OF U
%as explained in the constructive proof of the reference, loop over
%different target elements of U_current, the current version of U the
%target element is U_current(p,q)
for p=M:-1:2
    for q=p-1:-1:1
        %determine the current 2x2 submatrix Upq to be used determine the
        %current theta and phi values
        phi=atan2(imag(U_current(p,q))*real(U_current(p,p))-real(U_current(p,q))*imag(U_current(p,p)),imag(U_current(p,q))*imag(U_current(p,p))+real(U_current(p,q))*real(U_current(p,p)));
        theta=atan2(-(real(U_current(p,q))*cos(phi)+imag(U_current(p,q))*sin(phi)),(real(U_current(p,p))));
        %save these values for later
        phis(p,q)=phi;
        thetas(p,q)=theta;
        %from phi and theta, put together the 2x2 transformation
        T_qq=exp(-1i*phi)*cos(theta);
        T_qp=-exp(-1i*phi)*sin(theta);
        T_pq=sin(theta);
        T_pp=cos(theta);
        %construct Upq, the 2x2 transformation in the context of all the
        %modes
        Tpq=eye(M);
        Tpq(p,p)=T_pp;
        Tpq(p,q)=T_pq;
        Tpq(q,p)=T_qp;
        Tpq(q,q)=T_qq;
        %multiply from the right to update U_current
        U_current=U_current*Tpq;
    end
end

%correct the phases of the diagonal elements to make them equal to 1
phi2=NaN*ones(M,1);
for k=1:M
    phi2(k)=phase(U_current(k,k));
end
U_current=U_current*diag(exp(-1i*phi2));

%check that it worked (final U_current should be identity)
if norm(U_current-eye(M))>1E-8
    error('ReckDecomposition failed')
end

%% FORWARDS CONSTRUCTION OF U

%make sequence of simple transformation matrices
%do MZIs or simplifications thereof
for p=M:-1:2
    for q=p-1:-1:1
        Transformation_matrix_stack(:,:,m)=PhaseShifter(q,M,phis(p,q));
        m=m+1;
        Transformation_matrix_stack(:,:,m)=BeamSplitter(q,p,M,thetas(p,q));
        m=m+1;
    end
end

%and lastly the phases
for k=1:M
    if abs(phi2(k))>1E-8
        phase_incorporated='no';
        %can theta be incorporated in a previous phase?
        n=1;
        while (isdiag(Transformation_matrix_stack(:,:,max([m-n,1]))))&&(n<m)
            if (Transformation_matrix_stack(k,k,m-n)~=1) && (strcmp(phase_incorporated,'no'))
                Transformation_matrix_stack(k,k,m-n)=Transformation_matrix_stack(k,k,m-n)*exp(-1i*phi2(k));
                phase_incorporated='yes';
            end
            n=n+1;
        end
        if strcmp(phase_incorporated,'no')
            Transformation_matrix_stack(:,:,m)=PhaseShifter(k,M,phi2(k));
            m=m+1;
        end
    end
end

%deleting any identity operations to simplify and tidy up the result
deletenumber=[];
for n=1:size(Transformation_matrix_stack,3)
    if norm(Transformation_matrix_stack(:,:,n)-eye(M))<1E-8
        deletenumber(end+1,1)=n;
    end
end
if isempty(deletenumber)==0
    for m=size(deletenumber,1):-1:1
        Transformation_matrix_stack(:,:,deletenumber(m))=[];
    end
end

%multiply the final matrices to test that the construction works properly
U_test=eye(M);
for k=1:size(Transformation_matrix_stack,3)
    U_test=Transformation_matrix_stack(:,:,k)*U_test;
end
if norm(U_test-U)>1E-8
    error('something went wrong with the building blocks')
end
end

