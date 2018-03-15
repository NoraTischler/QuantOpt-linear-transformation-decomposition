function [ Transformation_matrix_stack ] = ClementsDecomposition( U,varargin )
%CLEMENTSDECOMPOSITION_VARBS provides the decomposition of a unitary matrix (U) into
%a sequence of basic 2x2 actions as outlined in the Optica paper
%http://dx.doi.org/10.1364/OPTICA.3.001460
%   The basic 2x2 actions are variable beam splitters (BS) and phase shifters
%   (PS). The output is given in terms of matrices that
%   describe the actions, in the order to be applied according to the third
%   dimension of the 3D array.

%% PRELIMINARIES AND INITIALISATION
% Check unitarity
tol=1E-8;
M=size(U,1);
if norm(U*U'-eye(M))+norm(U'*U-eye(M))>tol
    error('Input into ClementsDecomposition is not unitary')
end

%initialise variables
U_current=U;
T_left=eye(M);
phisleft=[];
thetasleft=[];
modes1=[];
modes2=[];
Transformation_matrix_stack_left=[];
m1=1;
m2=1;
m3=1;
m4=1;
phi2=NaN*ones(M,1);


%Default choice for the type of allowed beam splitters, if not specified
if nargin==1
    BStype='variable';
else 
    BStype=varargin{1};
end

%Default choice for the simplification if not specified
if nargin==1
    simplification='on';
else 
    simplification=varargin{2};
end
%% DECOMPOSITION OF U

%do a loop over target elements to null
for k1=1:M-1
    %nulling with T or Tinverse?
    if round(k1/2)==k1/2
        type='even'; %downwards nulling, multiplication from the left with T
    else
        type='odd'; %upwards nulling, multiplication from the right with Tinverse
    end
    
    switch type
        case 'even'
            TEr=M-k1;
            for TEc=1:k1
                %determine current target element
                TEr=TEr+1;
                %determine current transformation modes
                mode1=TEr-1;
                mode2=TEr;
                %find theta and phi values for the current target element
                phi=atan2(imag(U_current(mode1,TEc))*real(U_current(mode2,TEc))-real(U_current(mode1,TEc))*imag(U_current(mode2,TEc)),-imag(U_current(mode1,TEc))*imag(U_current(mode2,TEc))-real(U_current(mode1,TEc))*real(U_current(mode2,TEc)));
                theta=atan2(real(U_current(mode2,TEc)),cos(phi)*real(U_current(mode1,TEc))-sin(phi)*imag(U_current(mode1,TEc)));
                %record theta and phi values for the current target element
                phisleft(1,m1)=phi;
                thetasleft(1,m1)=theta;
                modes1(1,m1)=mode1;
                modes2(1,m1)=mode2;
                
                %from phi and theta, put together the 2x2 transformation
                T=eye(M);
                T(mode1,mode1)=exp(1i*phi)*cos(theta);
                T(mode2,mode1)=-exp(1i*phi)*sin(theta);
                T(mode1,mode2)=sin(theta);
                T(mode2,mode2)=cos(theta);
                %update U_current
                U_current=T*U_current;
                %update T_left
                T_left=T_left*T';
                m1=m1+1;
            end
        case 'odd'
            TEc=k1+1;
            for TEr=M:-1:M+1-k1
                %determine current target element
                TEc=TEc-1;
                %determine current transformation modes
                mode1=TEc;
                mode2=TEc+1;
                %find theta and phi values for the current target element
                phi=atan2(imag(U_current(TEr,mode1))*real(U_current(TEr,mode2))-real(U_current(TEr,mode1))*imag(U_current(TEr,mode2)),imag(U_current(TEr,mode1))*imag(U_current(TEr,mode2))+real(U_current(TEr,mode1))*real(U_current(TEr,mode2)));
                theta=atan2(-(real(U_current(TEr,mode1))*cos(phi)+imag(U_current(TEr,mode1))*sin(phi)),(real(U_current(TEr,mode2))));
                %record theta and phi values for the current target element
                phisright(1,m3)=phi;
                thetasright(1,m3)=theta;
                %from phi and theta, put together the 2x2 transformation
                Tinv=eye(M);
                Tinv(mode1,mode1)=exp(-1i*phi)*cos(theta);
                Tinv(mode1,mode2)=-exp(-1i*phi)*sin(theta);
                Tinv(mode2,mode1)=sin(theta);
                Tinv(mode2,mode2)=cos(theta);
                %update U_current
                U_current=U_current*Tinv;
                %update Transformation_matrix_stack_right
                switch BStype
                    case 'variable'
                        Transformation_matrix_stack_right(:,:,m4)=PhaseShifter(mode1,M,phisright(1,m3));
                        m4=m4+1;
                        Transformation_matrix_stack_right(:,:,m4)=BeamSplitter(mode1,mode2,M,thetasright(1,m3));
                        m4=m4+1;
                        m3=m3+1;
                    case 'balanced'
                        Transformation_matrix_stack_right(:,:,m4)=PhaseShifter(mode1,M,phisright(1,m3)+pi/2);
                        m4=m4+1;
                        Transformation_matrix_stack_right(:,:,m4)=BeamSplitter(mode1,mode2,M,pi/4);
                        m4=m4+1;
                        Transformation_matrix_stack_right(:,:,m4)=PhaseShifter(mode1,M,thetasright(1,m3)+pi);
                        m4=m4+1;
                        Transformation_matrix_stack_right(:,:,m4)=PhaseShifter(mode2,M,-thetasright(1,m3));
                        m4=m4+1;
                        Transformation_matrix_stack_right(:,:,m4)=BeamSplitter(mode1,mode2,M,pi/4);
                        m4=m4+1;
                        Transformation_matrix_stack_right(:,:,m4)=PhaseShifter(mode1,M,pi/2);
                        m4=m4+1;
                        m3=m3+1;
                end
            end
    end
end

phisleft=flip(phisleft,2);
thetasleft=flip(thetasleft,2);
modes1=flip(modes1,2);
modes2=flip(modes2,2);


%find D
D=U_current;
%calculate Dprime
Dprime=D;
for k=1:size(phisleft,2)
    %find commuted diagonal
    theta=thetasleft(1,k);
    phi=phisleft(1,k);
    psi1=phase(Dprime(modes1(1,k),modes1(1,k)));
    psi2=phase(Dprime(modes2(1,k),modes2(1,k)));
    thetaprime=-theta;
    phiprime=psi1-psi2;
    psi1prime=psi2-phi;
    psi2prime=psi2;
    phisleft(1,k)=phiprime;
    thetasleft(1,k)=thetaprime;
    Dprime(modes1(1,k),modes1(1,k))=exp(1i*psi1prime);
    Dprime(modes2(1,k),modes2(1,k))=exp(1i*psi2prime);
end

for k=1:M
    phi2(k)=phase(Dprime(k,k));
end
for k=1:M
    Transformation_matrix_stack_Dprime(:,:,k)=PhaseShifter(k,M,phi2(k));
    k=k+1;
end

    switch BStype
        case 'variable'
            for m1=1:size(phisleft,2)
                Transformation_matrix_stack_left(:,:,m2)=PhaseShifter(modes1(1,m1),M,phisleft(1,m1));
                m2=m2+1;
                Transformation_matrix_stack_left(:,:,m2)=BeamSplitter(modes1(1,m1),modes2(1,m1),M,thetasleft(1,m1));
                m2=m2+1;
            end
            
        case 'balanced'
            for m1=1:size(phisleft,2)
                Transformation_matrix_stack_left(:,:,m2)=PhaseShifter(modes1(1,m1),M,phisleft(1,m1)+pi/2);
                m2=m2+1;
                Transformation_matrix_stack_left(:,:,m2)=BeamSplitter(modes1(1,m1),modes2(1,m1),M,pi/4);
                m2=m2+1;
                Transformation_matrix_stack_left(:,:,m2)=PhaseShifter(modes1(1,m1),M,thetasleft(1,m1)+pi);
                m2=m2+1;
                Transformation_matrix_stack_left(:,:,m2)=PhaseShifter(modes2(1,m1),M,-thetasleft(1,m1));
                m2=m2+1;
                Transformation_matrix_stack_left(:,:,m2)=BeamSplitter(modes1(1,m1),modes2(1,m1),M,pi/4);
                m2=m2+1;
                Transformation_matrix_stack_left(:,:,m2)=PhaseShifter(modes1(1,m1),M,pi/2);
                m2=m2+1;
            end
    end
    

%% FORWARDS CONSTRUCTION OF U
Transformation_matrix_stack=cat(3,Transformation_matrix_stack_right,Transformation_matrix_stack_left,Transformation_matrix_stack_Dprime);

%simplifications if possible
switch simplification
    case 'on'
        Transformation_matrix_stack = Simplifications( Transformation_matrix_stack , M , 0 );
end

%multiply the final matrices to test that the construction works properly
U_test=eye(M);
for k=1:size(Transformation_matrix_stack,3)
    U_test=Transformation_matrix_stack(:,:,k)*U_test;
end
if norm(U_test-U)>tol
    error('Something went wrong with the building blocks')
end
end

