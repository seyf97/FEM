%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   In this function we update the ENL matrix based on the prescribed displacements.  %%%%%
%%%%   The columns number 9 and 10 correspond to the x- and y-component of the           %%%%%
%%%%   displacement at each node, respectively.                                          %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

function [ Up ]= assemble_displacements(ENL, NL, DOC)

NoN = size(NL,1);
Up = zeros(DOC,1);
PD = size(NL,2);

for i=1:NoN
    dummy = ENL(i,2*PD+1:3*PD); %temp
    
    if dummy(1)< 0
        Up(abs(dummy(1))) = ENL(i,4*PD+1);
    end
    
    if dummy(2)< 0
        Up(abs(dummy(2))) = ENL(i,4*PD+2);
    end
        
    
end





end
