function [ PS ] = PhaseShifter( mode , number_of_modes , phi )
%PHASESHIFTER creates the unitary matrix of a phase shift being implemented
%on the one specified mode, while the other modes are left alone
%   Inputs: 
%   mode, the mode on which the phase shifter is implemented
%   number_of_modes, the total number of modes
%   phi, the transformation parameter value specifying the phase shift
%   Output:
%   PS, the unitary matrix of the phase shifter transformation

PS=eye(number_of_modes);
PS(mode,mode)=exp(1i*phi);

end

