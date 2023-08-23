%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%   In this function we assign the boundary conditions required in the problem    %%%%%
%%%%   and we construct the ENL matrix based on the boundary conditions.             %%%%%
%%%%   Additionally we determine the degrees of freedom and degrees of constratint.  %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ ENL , DOF , DOC ] = assign_BCs( NL , Fgiven , Ugiven , BC )

NoN = size(NL,1); %1 is rows, 2 is columns
PD = size(NL,2); %Problem Dimension

%COORDINATES (1STBIG COLUMN)
ENL = zeros(NoN,6*PD);

ENL(:,1:PD) = NL;

%Initializing the BC, Disp, Force in ENL matrix

%BC info (2ND BIG COLUMN)
ENL(:,1*PD+1:2*PD ) = BC;

%Disp info (5TH BIG COLUMN)
ENL(:,4*PD+1:5*PD) = Ugiven;

%Force info (6TH BIG COLUMN)
ENL(:,5*PD+1:6*PD) = Fgiven;


%Initializing DoFs and DoCs  
DOF = 0;
DOC = 0;



%Forming the Temp Degree (3RD COLUMN)
Temp = zeros(NoN,PD);
for i = 1:NoN
    
    dummy = BC(i,:);
    
    %X Element
    if dummy(1) == -1
        DOC = DOC + dummy(1);
        Temp(i,1) = DOC;
    else
        DOF = DOF + dummy(1);
        Temp(i,1) = DOF;
    end
    
    %Y Element
    if dummy(2) == -1
        DOC = DOC + dummy(2);
        Temp(i,2) = DOC;
    else
        DOF = DOF + dummy(2);
        Temp(i,2) = DOF;
    end
end

%Now we have the Dof and Doc numbers. We take the abs value of the DOC
DOC = abs(DOC);

%Finding the global degrees. + Temp elements stay the same (4TH COLUMN):
Global = Temp;

for i = 1:NoN
    
    dummy = Temp(i,:);
    
    %X Element
    if dummy(1) < 0
        Global(i,1) = abs(dummy(1)) + DOF;
    end
    
    %Y Element
    if dummy(2) < 0
        Global(i,2) = abs(dummy(2)) + DOF;
    end
end

ENL(:,2*PD+1:3*PD) = Temp;
ENL(:,3*PD+1:4*PD) = Global;



end