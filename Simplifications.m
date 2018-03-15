function [ S_stack_total_simplified ] = Simplifications( S_stack_total , n , nA )
%SIMPLIFICATIONS takes the sequence of transformations in the form or
%S_stack_total as an input, and performs a number of simplifications if
%possible.
%   Mach-Zehnder interferometers are replaced by simpler elements, swap
%   operations are permuted to the end, phase shifters are combined, and
%   any identity operations are deleted.
%   Inputs:
%   S_stack_total, the sequence of elementary transformations
%   n, the number of nominal modes
%   nA, the number of ancilla modes
%   Outputs: 
%   S_stack_total_simplified, the simplified version of the sequence of
%   transformations

tol=1E-8;
if size(S_stack_total,1)==2*(n+nA)
    MatrixType='S';
elseif size(S_stack_total,1)==n+nA
    MatrixType='U';
    S_stack_temp=repmat(eye(2*(n+nA)),[1, 1, size(S_stack_total,3)]);
    S_stack_temp(1:size(S_stack_total,1),1:size(S_stack_total,1),:)=S_stack_total;
    S_stack_total=S_stack_temp;
else
    error('The number of rows and columns of the first input into Simplifications must be either equal to n+nA, or equal to 2*(n+nA).')
end
%simplify any Mach-Zehnder interferometers where the two beam splitters
%are not needed
for operation_number=1:size(S_stack_total,3)-3
    simp=0;
    [mode1,~]=find(abs(S_stack_total(1:n+nA,1:n+nA,operation_number)-(1/sqrt(2)))<tol,1,'first');
    [mode2,~]=find(abs(S_stack_total(1:n+nA,1:n+nA,operation_number)-(1/sqrt(2)))<tol,1,'last');
    [mode1b,~]=find(abs(S_stack_total(1:n+nA,1:n+nA,operation_number+3)-(1/sqrt(2)))<tol,1,'first');
    [mode2b,~]=find(abs(S_stack_total(1:n+nA,1:n+nA,operation_number+3)-(1/sqrt(2)))<tol,1,'last');
    if (isempty(mode1)==0&&isempty(mode1b)==0)
        if (mode1b==mode1 && mode2b==mode2) %identify locations of MZIs
            
            MZIop=S_stack_total(:,:,operation_number+3)*S_stack_total(:,:,operation_number+2)*S_stack_total(:,:,operation_number+1)*S_stack_total(:,:,operation_number);
            if norm(abs(MZIop([mode1 mode2],[mode1 mode2]))-eye(2))<tol %check if it falls into one of the categories of identity, swap, or single beam splitter
                trans=eye(2);
                resid=MZIop([mode1 mode2],[mode1 mode2])/trans; %calculate any necessary phase shifts or swaps
                simp=1;
            elseif norm(abs(MZIop([mode1 mode2],[mode1 mode2]))-[0,1;1,0])<tol
                trans=[0,1;1,0];
                resid=MZIop([mode1 mode2],[mode1 mode2])/trans;
                simp=1;
            elseif norm(abs(MZIop([mode1 mode2],[mode1 mode2]))-(1/sqrt(2)*[1,1;1,1]))<tol
                trans=1/sqrt(2)*[1,1;-1,1];
                resphases=angle(MZIop([mode1 mode2],[mode1 mode2])./trans);
                if abs(resphases(2,2)-(resphases(1,2)+resphases(2,1)-resphases(1,1)))<tol
                    phase2pre=-resphases(1,1)+resphases(1,2);
                    phase1post=resphases(1,1);
                    phase2post=resphases(2,1);
                    UM=eye(n+nA);
                    UM(mode2,mode2)=exp(1i*phase2pre);
                    S_stack_total(:,:,operation_number)=EnlargeUnitary(UM,0);
                    UM([mode1 mode2],[mode1 mode2])=trans;
                    S_stack_total(:,:,operation_number+1)=EnlargeUnitary(UM,0);
                    UM=eye(n+nA);
                    UM(mode1,mode1)=exp(1i*phase1post);
                    S_stack_total(:,:,operation_number+2)=EnlargeUnitary(UM,0);
                    UM=eye(n+nA);
                    UM(mode2,mode2)=exp(1i*phase2post);
                    S_stack_total(:,:,operation_number+3)=EnlargeUnitary(UM,0);
                    simp=0;
                end
            end
            if simp==1
                if isdiag(round(resid,8))
                    simp=1;
                elseif norm(abs(resid)-[0,1;1,0])<tol
                    simp=2;
                else
                    simp=0;
                end
            end
            if simp==1 %substitute simplified elements in place of the MZI
                S_stack_total(:,:,operation_number:operation_number+3)=repmat(eye(size(S_stack_total,1)),[1 1 4]);
                UM=eye(n+nA);
                UM([mode1 mode2],[mode1 mode2])=trans;
                S_stack_total(:,:,operation_number)=EnlargeUnitary(UM,0);
                UM([mode1 mode2],[mode1 mode2])=diag([resid(1,1),abs(resid(2,2))]);
                S_stack_total(:,:,operation_number+1)=EnlargeUnitary(UM,0);
                UM([mode1 mode2],[mode1 mode2])=diag([abs(resid(1,1)), resid(2,2)]);
                S_stack_total(:,:,operation_number+3)=EnlargeUnitary(UM,0);
            elseif simp==2
                S_stack_total(:,:,operation_number:operation_number+3)=repmat(eye(size(S_stack_total,1)),[1 1 4]);
                UM([mode1 mode2],[mode1 mode2])=trans;
                S_stack_total(:,:,operation_number)=EnlargeUnitary(UM,0);
                UM([mode1 mode2],[mode1 mode2])=resid;
                S_stack_total(:,:,operation_number+1)=EnlargeUnitary(UM,0);
                if abs((sum(sum(resid))+2) < tol)
                    UM([mode1 mode2],[mode1 mode2])=[0,1;-1,0];
                    S_stack_total(:,:,operation_number+1)=EnlargeUnitary(UM,0);
                    UM([mode1 mode2],[mode1 mode2])=[-1,0;0,1];
                    S_stack_total(:,:,operation_number+2)=EnlargeUnitary(UM,0);
                end
            end
        end
    end
end
%commute swap operations (permutations of modes) to the end of the setup
for operation_number=1:size(S_stack_total,3)
    swap_modes=find(abs(diag(S_stack_total(:,:,operation_number)))<tol,2,'first');%find swap operation
    swap_modes2=find(abs(diag(S_stack_total(:,:,operation_number)))<tol,2,'last');
    if isempty(swap_modes)==0
        perm_matrix=S_stack_total(:,:,operation_number);
        [resid_phase,neg_el_col]=min([S_stack_total(swap_modes(2),swap_modes(1),operation_number),S_stack_total(swap_modes(1),swap_modes(2),operation_number)]);
        phase_mode=swap_modes(neg_el_col);
        S_stack_total(:,:,operation_number)=eye(size(S_stack_total,1));
        if (resid_phase+1)<tol %check if there was a phase of pi, or if it was a simple swap
            S_stack_total(phase_mode,phase_mode,operation_number)=-1; %implement the pi phase shift in place of the swap
            S_stack_total(phase_mode+n+nA,phase_mode+n+nA,operation_number)=-1;
        end
        perm_vec=linspace(1,size(S_stack_total,1),size(S_stack_total,1));
        perm_vec(swap_modes) = perm_vec(flip(swap_modes,1));
        perm_vec2=linspace(1,size(S_stack_total,1),size(S_stack_total,1));
      %  perm_vec2(swap_modes2) = perm_vec(flip(swap_modes2,1));%make permutation vector
      if swap_modes~=swap_modes2
        perm_vec2(swap_modes2) = perm_vec2(flip(swap_modes2,1));%make permutation vector
      end
        for k=operation_number+1:size(S_stack_total,3)
            S_stack_total(:,:,k) = S_stack_total(:,perm_vec,k);
            S_stack_total(:,:,k) = S_stack_total(perm_vec,:,k);
            S_stack_total(:,:,k) = S_stack_total(:,perm_vec2,k);
            S_stack_total(:,:,k) = S_stack_total(perm_vec2,:,k);%propagate the permutation through to the end
        end
  %      S_stack_total(:,:,end+1)=perm_matrix;%place the permutation at the end
 S_stack_total(:,:,end+1)=abs(perm_matrix);%place the permutation at the end
    end
end

%identify double phase shifts from BS operation 
for operation_number=1:2*size(S_stack_total,3)
    if operation_number<=size(S_stack_total,3)
        phases=find((abs(diag(S_stack_total(1:n+nA,1:n+nA,operation_number))-1)>tol) & abs(abs(diag(S_stack_total(1:n+nA,1:n+nA,operation_number)))-1)<tol);
        if length(phases)>1
            S_stack_total(:,:,operation_number+1:size(S_stack_total,3)+1)=S_stack_total(:,:,operation_number:end);
            S_stack_total(phases(2),phases(2),operation_number)=1;
            S_stack_total(phases(1),phases(1),operation_number+1)=1;
            S_stack_total(phases(2)+n+nA,phases(2)+n+nA,operation_number)=1;
            S_stack_total(phases(1)+n+nA,phases(1)+n+nA,operation_number+1)=1;
        end
    end
end

%combine phase shifters if possible
for mode_number=1:n+nA
    indBS=find(abs(abs(S_stack_total(mode_number,mode_number,:))-1)>tol);
    indPS=find((abs(S_stack_total(mode_number,mode_number,:)-1)>tol) & abs(abs(S_stack_total(mode_number,mode_number,:))-1)<tol);%make vectors of phase shift locations and of beam splitter
    %locations
    indgroups = discretize(indPS,[0; indBS; size(S_stack_total,3)+1]); %groups of phase shifts that can be combined
    for k1=length(indgroups):-1:1
        curr_group=indgroups(k1);
        PStocombine=find(indgroups==curr_group,2,'last');
        if length(PStocombine)>1
            %combine phases
            S_stack_total(:,:,indPS(PStocombine(1,1)))=S_stack_total(:,:,indPS(PStocombine(1,1))).*S_stack_total(:,:,indPS(PStocombine(2,1)));
            S_stack_total(:,:,indPS(PStocombine(2,1)))=[];
            indgroups(PStocombine(1,1))=[];
        end
    end
end
%delete any identity operations to simplify and tidy up the result
deletenumber=[];
for k=1:size(S_stack_total,3)
    if norm(S_stack_total(:,:,k)-eye(size(S_stack_total,1)))<tol
        deletenumber(end+1,1)=k;
    end
end
if isempty(deletenumber)==0
    for m=size(deletenumber,1):-1:1
        S_stack_total(:,:,deletenumber(m))=[];
    end
end

switch MatrixType
    case 'S'
        S_stack_total_simplified=S_stack_total;
    case 'U'
        S_stack_total_simplified=S_stack_total(1:n+nA,1:n+nA,:);
end
end

