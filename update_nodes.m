%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   In this function we update the ENL matrix based on the solution of      %%%%
%%%%   our problem. The unknown displacement and forces that have been sought  %%%%  
%%%%   are put in the correct position.                                        %%%%       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  ENL  = update_nodes( ENL , U, NL)

NoN = size(NL,1);
PD = size(NL,2);

for i=1:NoN
    
    dummy = ENL(i,3*PD+1:4*PD);
    
    %Putting the displacements
    ENL(i,4*PD+1) = U(dummy(1));
    ENL(i,5*PD)=    U(dummy(2));
    
%     %Putting the forces
%     ENL(i,5*PD+1) = F(dummy(1));
%     ENL(i,6*PD)=    F(dummy(2));


end