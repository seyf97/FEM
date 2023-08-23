%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   In this function we update the ENL matrix based on the prescribed forces.  %%%%%
%%%%   The columns number 11 and 12 correspond to the x- and y-component of the   %%%%%
%%%%   force vectors at each node, respectively.                                  %%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Fp = assemble_forces( ENL , NL, DOF)

NoN = size(NL,1);

Fp = zeros(DOF,1);

PD = size(NL,2);

for i=1:NoN
    dummy = ENL(i,2*PD+1:3*PD); %temp
    
    if dummy(1)> 0
        Fp(dummy(1)) = ENL(i,5*PD+1);
    end
    
    if dummy(2)>0
        Fp(dummy(2)) = ENL(i,5*PD+2);
    end
        
    
    
end




end