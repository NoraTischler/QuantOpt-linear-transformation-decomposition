function [ Swap ] = Swapper( mode1, mode2,number_of_modes )
%SWAPPER creates the unitary matrix where a swap is implemented on the
%two specified modes, and the other modes are left alone

Swap=eye(number_of_modes);
Swap(mode1,mode1)=0;
Swap(mode1,mode2)=1;
Swap(mode2,mode1)=1;
Swap(mode2,mode2)=0;


end

