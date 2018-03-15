function [ S_total, S_stack_total,Transformation_summary ] = GeneralTransformationSolver( T )
%GENERALTRANSFORMATIONSOLVER Takes the general linear transformation T, and
%computes a corresponding full transformation (S matrix), the simple S
%matrices it can be decomposed into, and the corresponding sequence of
%simple transformations. The allowed building blocks are phase shifters,
%beam splitters (either variable or balanced, as specified in the options
%section), and parametric amplifiers.

%   Input: 
%   T, a complex matrix of any size, specifying the linear transformation
%   of interest

%   Outputs: S_total: A full bosonic transformation including ancilla modes
%   (a complex matrix including T as its first block) S_stack_total: A 3D
%   array providing the bosonic transformations of the step-by-step setup
%   (2D transformation matrices, with the step indexed by the third
%   dimension of the array) Transformation_summary: A table detailing the
%   sequence of elemental transformations comprising the setup. The table
%   contains the type of transformation, the value of the parameter*pi
%   (rad), and the modes between which the transformation takes place. The
%   corresponding 1x1 or 2x2 transformations for a parameter value of x
%   are: phase shifter = [exp(i*x)] beam splitter =
%   [cos(x),sin(x);-sin(x),cos(x)] swap = [0,1;1,0] identity = [1,0;0,1]
%   parametric amplifier =
%   [a1_out;a2_out;a1_out^dagger;a2_out^dagger]...
%   ...=[cosh(x),0,0,sinh(x);0,cosh(x),sinh(x),0;...
%   ...0,sinh(x),cosh(x),0;sinh(x),0,0,cosh(x)]...
%   ...*[a1_in;a2_in;a1_in^dagger;a2_in^dagger]

%   Options: The type of beam splitter that is allowed (either variable
%   splitting ratio or balanced), whether the setup should be simplified
%   where possible or kept according to the original decomposition, and
%   which of the two available decomposition methods for unitary matrices
%   should be used

%% OPTIONS
%pick one option for each pair below

BStype='variable'; %variable beam splitting ratio allowed
%BStype='balanced'; %only 50-50 beam splitting ratio allowed

%simplification='off'; %suited to integrated optics implementation with a univeral circuit containing tunable phase shifters
simplification='on'; %suited to bulk optics implementation, where the circuit should be as simple as possible

unitary_decomposition_method='Reck'; % the unitary decomposition follows the method from https://doi.org/10.1103/PhysRevLett.73.58
%unitary_decomposition_method='Clements'; % the unitary decomposition follows the method from http://dx.doi.org/10.1364/OPTICA.3.001460

%% DECOMPOSITION

%preliminaries and initialization
tol=1E-6;
n=max(size(T));
nA=0;
U=eye(n);
D=eye(n);
Wt=eye(n);

%do singular value decomposition (SVD)
[U(1:size(T,1),1:size(T,1)),D(1:size(T,1),1:size(T,2)),Wt(1:size(T,2),1:size(T,2))]=svd(T);
W=Wt';

%determine number of ancillas
for k=1:n
    if abs(D(k,k)-1)>tol
        nA=nA+1;
    end
end

%S matrices for U
switch unitary_decomposition_method
    case 'Reck'
        Transformation_matrix_stack_U = ReckDecomposition( U,BStype,simplification ); 
    case 'Clements'
        Transformation_matrix_stack_U = ClementsDecomposition( U,BStype,simplification ); 
end
S_stack_U=NaN*ones(2*(n+nA),2*(n+nA),size(Transformation_matrix_stack_U,3));
if size(Transformation_matrix_stack_U,3)>0
    for k=1:size(Transformation_matrix_stack_U,3)
        S_stack_U(:,:,k)=EnlargeUnitary( Transformation_matrix_stack_U(:,:,k),nA );
    end
end

%S matrices for V
switch unitary_decomposition_method
    case 'Reck'
        Transformation_matrix_stack_V = ReckDecomposition( W,BStype,simplification );
    case 'Clements'
        Transformation_matrix_stack_V = ClementsDecomposition( W,BStype,simplification );
end
S_stack_V=NaN*ones(2*(n+nA),2*(n+nA),size(Transformation_matrix_stack_V,3));
if size(Transformation_matrix_stack_V,3)>0
    for k=1:size(Transformation_matrix_stack_V,3)
        S_stack_V(:,:,k)=EnlargeUnitary( Transformation_matrix_stack_V(:,:,k),nA );
    end
end

%S matrices for D
k=1;
pos_ancilla=1;
S_stack_D=[];
for pos_el=1:n
    if abs(D(pos_el,pos_el)-1)>tol
        S_p_to_fill=S_p_from_diag( D, pos_el,pos_ancilla,nA,BStype );
        S_stack_D(:,:,k:k+size(S_p_to_fill,3)-1)= S_p_to_fill;
        k=size(S_stack_D,3)+1;
        pos_ancilla=pos_ancilla+1;
    end
end

%a step-by-step decomposition in order of transformation application
S_stack_total=cat(3,S_stack_V,S_stack_D,S_stack_U);

%simplifications of the setup if possible and desired
switch simplification
    case 'on'
    [S_stack_total] = Simplifications(S_stack_total,n,nA);
end

%multiply the final matrices to get the full overall transformation
S_total=eye(2*(n+nA));
for k=1:size(S_stack_total,3)
    S_total=S_stack_total(:,:,k)*S_total;
end

%% TESTING

Ssize=size(S_total,1);
%test quasiunitarity of S_total
G=MakeG( Ssize );
if norm(S_total*G*S_total'-G)>tol
    error('S_total is not quasiunitary')
end
%test that the construction works properly, such that T is the first block
%of S_total
if norm(S_total(1:size(T,1),1:size(T,2))-T)>tol
    error('Something went wrong with the building blocks')
end

%% GENERATING A SUMMARY OF THE SETUP

Transformation_type=NaN*ones(size(S_stack_total,3),1);
First_Mode=NaN*ones(size(S_stack_total,3),1);
Second_Mode=NaN*ones(size(S_stack_total,3),1);
Parameter_Value=NaN*ones(size(S_stack_total,3),1);

for k=1:size(S_stack_total,3)
    [mode_involved,~]=find(abs(diag(diag(S_stack_total(:,:,k)))-eye(size(S_stack_total,2)))>tol);
    if isempty(mode_involved)
        Transformation_type(k)=0;
    else
        First_Mode(k)=mode_involved(1);
        if length(mode_involved)>2
            Second_Mode(k)=mode_involved(2);
            if norm(abs(S_stack_total([mode_involved(1) mode_involved(2)],[mode_involved(1) mode_involved(2)],k)-[0,1;1,0]))<tol
                Transformation_type(k)=4;
            else
                if S_stack_total(mode_involved(1),mode_involved(1),k)<1
                    Transformation_type(k)=2;
                    Parameter_Value(k)=atan2(S_stack_total(mode_involved(1),mode_involved(2),k),S_stack_total(mode_involved(1),mode_involved(1),k))/pi;
                    
                elseif S_stack_total(mode_involved(1),mode_involved(1),k)>1
                    Transformation_type(k)=3;
                    Parameter_Value(k)=acosh(S_stack_total(mode_involved(1),mode_involved(1),k))/pi;
                end
            end
        else
            Transformation_type(k)=1;
            Parameter_Value(k)=angle(S_stack_total(mode_involved(1),mode_involved(1),k))/pi;
        end
    end
end
Transformation_type= num2cell(Transformation_type);
for k=1:size(Transformation_type,1)
    if Transformation_type{k}==0
        Transformation_type{k}='identity';
    elseif Transformation_type{k}==1
        Transformation_type{k}='phase shifter';
    elseif Transformation_type{k}==2
        Transformation_type{k}='beam splitter';
    elseif Transformation_type{k}==3
        Transformation_type{k}='parametric amplifier';
    elseif Transformation_type{k}==4
        Transformation_type{k}='mode swap';
    end
end
Transformation_type = categorical(Transformation_type);
Transformation_summary = table(Transformation_type,Parameter_Value,First_Mode,Second_Mode); %the parameter value is given in units of pi radians
end

