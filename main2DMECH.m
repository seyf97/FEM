clc
clear
close all

% PD : Problem Dimension
% NPE : Number of Nodes per Element

% NoN : Number of Nodes
% NoE : Number of Elements

% NL : Nodes List [ NoN X PD ]
% EL : Elements List [ NoE X NPE ]

% ENL : Extended Node List [ NoN X (6xPD) ]

% DOFs : Degrees of Freedom
% DOCs : Degrees of Constraint

% Up : Prescribed Displacements [ DOCs X 1 ]
% U : Unknown Displacements [ DOFs X 1 ]

% K_reduced : Reduced Stiffness [ DOFs X DOFs ]
% Fp : Prescribed Forces [ DOFs X 1 ]
% F : Right-hand Sided [ DOFs X 1 ]

% mag : Magnification Factor (FOR POST PROCESSING PART)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%   PRE-PROCESS.   %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format long

r = 0.2; %same for cube and rhombus
lengthy = 1;
l = lengthy;
lengthx = 1;
w = lengthx;
elType = 'D2TR3N';
inclType = 'C'; %(R)hombus or (S)quare as well as (C)ircle 
inclFilled = false; %is inside filled
GPE = 4;    
numofsegy = 15;
xmd = 15;
cgx = 0;
cgy = 0;

Eo = 80/3;
nuO = 1/3;

Ei = 10;
nuI = 0.3;

[EL,NL] = meshFinal(r,l,w,elType,numofsegy,xmd,cgx,cgy,inclType,inclFilled);


NPE = size(EL,2);

hold on
for i = 1:size(EL,1)

    Nodes = EL(i,:);
%     a = (NL(Nodes(1),:) + NL(Nodes(2),:) + NL(Nodes(3),:))./3;
%     s = sprintf('%d',i);
%     text(a(1),a(2),s)
% 
%     scatter(a(1),a(2),'k')
    patch(NL(Nodes,1),NL(Nodes,2),0,'FaceColor','none')
end

scatter(NL(:,1),NL(:,2),'k','filled')
axis equal
grid on

[ELnew,NLnew] = advancedElements(elType,NL,EL);


NoN = size(NL,1); %1 is rows, 2 is columns
PD = size(NL,2); %Problem Dimension


BC = ones(NoN,PD); %Initializing BC for all f = 0 BC's
%Initializing prescribred Force matrix
Fgiven = zeros(NoN,PD);

%Initializing prescribed Disp. matrix
Ugiven = zeros(NoN,PD);

%     PRESCRIBED DISPLACEMENTS    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%EDIT THIS PART FOR EVERY PROBLEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%Shear Test
% for i = 1:NoN
%     
%     if round(NL(i,2),5) == 0
%         BC(i,:) = [-1 -1];
%         Ugiven(i,:) = [0 0];
% %         scatter(NL(i,1),NL(i,2),'r')
%         
%     elseif round(NL(i,2),5) == l
%         BC(i,:) = [-1 -1];
%         Ugiven(i,:) = [0.1 0];
% %         scatter(NL(i,1),NL(i,2),'b')
% 
%     end
%     
%     
% end


% % %%%%%%%%%%%%%%%%%%Expansion Test%%%%%%%%%%%%%%%%%%
% for i = 1:NoN
%     
%     if round(NL(i,1),5) == 0 || round(NL(i,1),5) == l || round(NL(i,2),5) == 0 ...
%             || round(NL(i,2),5) == l
%         BC(i,:) = [-1 -1];
%         Ugiven(i,:) = [0.1*NL(i,1) 0.1*NL(i,2)];
%     end
%     
% end
BC(end,:) = [-1 -1];
Ugiven(end,:) = [0.1 0.1]; 
% % 
% % %%%%%%%%%%%%%%%%%%%%%% Extension Test
% % for i = 1:NoN
% %     
% %     if round(NL(i,1),5) == 0
% %         BC(i,:) = [-1 -1];
% %         Ugiven(i,:) = [0 0];
% %         
% %     elseif round(NL(i,1),5) == lengthy
% %         BC(i,:) = [-1 -1];
% %         Ugiven(i,:) = [0.2 0];
% %         
% %     end
% %     
% % end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%   PROCESS    %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[ENL,DOF,DOC] = assign_BCs(NL,Fgiven,Ugiven,BC);

K = GlobalStiffness(EL,NL,ENL,GPE,Ei,nuI,Eo,nuO,numofsegy,xmd,NPE);
% E_2 = 0;
% nu_2 = 0;
% E_1 = 1;
% nu_1 = 0.3;
% 
% K2 = assemble_stiffness(ENL,EL,NL,E_1,nu_1,E_2,nu_2);


Fp = assemble_forces( ENL , NL, DOF);

Up = assemble_displacements(ENL,NL,DOC);

K_11 = K(1:DOF,1:DOF);
K_12 = K(1:DOF,DOF+1:end);
K_21 = K(DOF+1:end,1:DOF);
K_22 = K(DOF+1:end,DOF+1:end);

G = Fp - K_12* Up; %G

Uu = K_11 \ G; %[A]^-1*G BU COMMANDI KULLAN
% 
Fu = K_21*Uu + K_22*Up;
% 
F = [Fp;Fu];
U = [Uu;Up];

ENL  = update_nodes( ENL , U, NL);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%   POST-PROCESS    %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mag = 1 --> Original
% Mag = 8000;
%8000
[NoN] = post_process2(NL,EL,ENL,Ei,nuI,Eo,nuO,numofsegy,xmd,GPE);

% [Stress, Strain, NewCoord,L0,Lf] = findStressStrain(ENL, E, NoN, EL, PD);

