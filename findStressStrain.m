function [Stress, Strain, NewCoord,L0,Lf] = findStressStrain(ENL, E, NoN, EL, PD)

%New Coordinates of the nodes
NewCoord = ENL(:,1:2) + ENL(:,9:10);
L0 = zeros(NoN,1);
Lf = zeros(NoN,1);

Stress = zeros(size(EL,1),3);
Strain = zeros(size(EL,1),3);

for i = 1:NoN
    
    Strain(i,1) = ENL(i,9); %xx
    Strain(i,2) = ENL(i,10); %yy
    Strain(i,3) = ENL()
    
end


    
    
    
end




    