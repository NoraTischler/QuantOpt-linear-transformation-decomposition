function [ BS ] = BeamSplitter( mode1 , mode2 , number_of_modes , theta )
%BEAMSPLITTER creates the unitary matrix of a variable beam splitter being
%implemented on the specified two modes, while the other modes are left
%alone
%   Inputs: 
%   mode1, the first mode of the beam splitter
%   mode2, the first mode of the beam splitter
%   number_of_modes, the total number of modes
%   theta, the transformation parameter value specifying the splitting
%   ratio
%   Output:
%   BS, the unitary matrix of the beam splitter transformation

BS=eye(number_of_modes);
BS(mode1,mode1)=cos(theta);
BS(mode1,mode2)=sin(theta);
BS(mode2,mode1)=-sin(theta);
BS(mode2,mode2)=cos(theta);

end

