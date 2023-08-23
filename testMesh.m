clc
clear
close all

r = 0.2;
lengthy = 1;
l = lengthy;
lengthx = 1;
w = lengthx;
elType = 'D2QU4N';
inclType = 'C'; %(R)hombus or (S)quare as well as (C)ircle 
inclFilled = true; %is inside filled
GPE = 3;    
numofsegy = 4;
xmd = 2;
cgx = 0;
cgy = 0;


[EL,NL] = meshFinal(r,l,w,elType,numofsegy,xmd,cgx,cgy,inclType,inclFilled);



hold on
axis equal
for i = 1:size(EL,1)

    Nodes = EL(i,:);
    a = (NL(Nodes(1),:) + NL(Nodes(2),:) + NL(Nodes(3),:))./3;
    s = sprintf('%d',i);
    text(a(1),a(2),s)

    scatter(a(1),a(2),'k')
    plot(NL([Nodes Nodes(1)],1),NL([Nodes Nodes(1)],2),'k')

end




function [EL,NL] = meshFinal(r,l,w,elType,numofsegy,xmd,cgx,cgy,inclType,inclFilled)

lengthx = w;
lengthy = l;
numofsegc = 4*(numofsegy);



             %%%%%%%%%%%%%%%%         1          %%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%  OUTER RECTANGLE   %%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%                    %%%%%%%%%%%%%%%%

xl = linspace(-l/2,l/2,numofsegy+1);
xw = linspace(-w/2,w/2,numofsegy+1);

yw = zeros(1,numofsegy+1);
yl = zeros(1,numofsegy+1);

xb = [xw' (yw'-l/2)];
xr = [(yl'+w/2) xl'];

xu = xb;
xu(:,2) = xu(:,2) + l;
xl = xr;
xl(:,1) = xl(:,1) - w;

x0 = [xb;xr(2:end,:);xu(end-1:-1:1,:);xl(end-1:-1:2,:)];
steps = numofsegy/2;
dummy = x0;

if  inclType == 'R'
    %Shifting the first point
    x0(1:end-steps,:) = x0(1+steps:end,:);
    x0(end-steps+1:end,:) = dummy(1:steps,:);
end
x0 = x0 + [lengthx/2 lengthy/2]; %moving to 00




             %%%%%%%%%%%%%%%%         2          %%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%     INCLUSION      %%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%                    %%%%%%%%%%%%%%%%




if inclType == 'C' %INCLUSION: C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xc = r.*cosd(linspace(225,360+225,numofsegc+1));
    yc = r.*sind(linspace(225,360+225,numofsegc+1));
    %shifting the CG of the Circle
    xc = xc + cgx ;
    yc = yc + cgy;
    
    xc(end) = [];
    yc(end) = [];
    xc = [xc' yc'];
    
    xc = xc + [lengthx/2 lengthy/2]; %moving to center
    xi = xc;

    if inclFilled
        
        xl = linspace(-r/(3*sqrt(2)),r/(3*sqrt(2)),numofsegy+1);
        xw = linspace(-r/(3*sqrt(2)),r/(3*sqrt(2)),numofsegy+1);
        
        yw = zeros(1,numofsegy+1);
        yl = zeros(1,numofsegy+1);
        
        xb = [xw' (yw'-r/(3*sqrt(2)))];
        xr = [(yl'+r/(3*sqrt(2))) xl'];
        
        xu = xb;
        xu(:,2) = xu(:,2) + 2*r/(3*sqrt(2));
        xl = xr;
        xl(:,1) = xl(:,1) - 2*r/(3*sqrt(2));
        
        x01 = [xb;xr(2:end,:);xu(end-1:-1:1,:);xl(end-1:-1:2,:)];
        x01 = x01 + [lengthx/2 lengthy/2]; %moving to center
        xi2 = x01;
        
        %Corners
        x1 = xi2(1,:);
        x2 = xi2(1*numofsegy+1,:);
        x3 = xi2(2*numofsegy+1,:);
        x4 = xi2(3*numofsegy+1,:);
        
    end
    
    
    
    
    
    
elseif inclType == 'R' %INCLUSION RHOMBUS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xl = linspace(-r/2,r/2,numofsegy+1);
    xw = linspace(-r/2,r/2,numofsegy+1);
    
    yw = zeros(1,numofsegy+1);
    yl = zeros(1,numofsegy+1);
    
    xb = [xw' (yw'-r/2)];
    xr = [(yl'+r/2) xl'];
    
    xu = xb;
    xu(:,2) = xu(:,2) + r;
    xl = xr;
    xl(:,1) = xl(:,1) - r;
    x0r = [xb;xr(2:end,:);xu(end-1:-1:1,:);xl(end-1:-1:2,:)];
    
    
    %rotating to obtain the rhombus
    newAngles = pi/4 + atan2(x0r(:,2),x0r(:,1));
    lengths = sqrt(x0r(:,1).^2 + x0r(:,2).^2);
    x0r = [lengths.*cos(newAngles) lengths.*sin(newAngles)];
    
    x0r = x0r + [lengthx/2 lengthy/2]; %moving to center
    
    steps = numofsegy/2;
    dummy = x0r;
    
    %Corners
    x1 = x0r(1,:);
    x2 = x0r(1*numofsegy+1,:);
    x3 = x0r(2*numofsegy+1,:);
    x4 = x0r(3*numofsegy+1,:);
    
    xi = x0r;
    
    
elseif inclType == 'S' %INCLUSION SQUARE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    xl = linspace(-r/2,r/2,numofsegy+1);
    xw = linspace(-r/2,r/2,numofsegy+1);
    
    yw = zeros(1,numofsegy+1);
    yl = zeros(1,numofsegy+1);
    
    xb = [xw' (yw'-r/2)];
    xr = [(yl'+r/2) xl'];
    
    xu = xb;
    xu(:,2) = xu(:,2) + r;
    xl = xr;
    xl(:,1) = xl(:,1) - r;
    xi = [xb;xr(2:end,:);xu(end-1:-1:1,:);xl(end-1:-1:2,:)];
    
    xi = xi + [lengthx/2 lengthy/2]; %moving to center

    %Corners
    x1 = xi(1,:);
    x2 = xi(1*numofsegy+1,:);
    x3 = xi(2*numofsegy+1,:);
    x4 = xi(3*numofsegy+1,:);
end



             %%%%%%%%%%%%%%%%         3          %%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%      MESHING       %%%%%%%%%%%%%%%%
             %%%%%%%%%%%%%%%%                    %%%%%%%%%%%%%%%%


%NL1 AND EL1
%AREA BTW INCLUSION AND OUTER RECTANGLE. NL1
NL1 = zeros((xmd+2)*numofsegc,2);
cnt = 1;
for i = 1:numofsegc
    
    dummy = lineGen(x0(i,:),xi(i,:),xmd+1);
    dummy(end,:) = [];
    NL1(cnt:cnt+size(dummy,1)-1,:) = dummy;
    cnt = cnt + size(dummy,1);
    
end
NL1(cnt:end,:) = xi;
EL1 = middleMesh(elType,NL1,numofsegy,xmd);





if inclType == 'C' && inclFilled
    
    %AREA BTW INCLUSION AND OUTER Circle. NL11
    NL11 = zeros((numofsegc)*numofsegc,2); %RANDOM
    cnt = 1;
    for i = 1:numofsegc
        dummy = lineGen(xi(i,:),x01(i,:),xmd+1);
        dummy(end,:) = [];
        dummy(1,:) = [];
        NL11(cnt:cnt+size(dummy,1)-1,:) = dummy;
        cnt = cnt + size(dummy,1);
    end
    NL11(cnt:cnt+size(xi2,1)-1,:) = xi2;
    NL11 = [xi;NL11];
    cnt = cnt + size(xi2,1) + size(xi,1);
    NL11(cnt:end,:) = [];
    EL11 = middleMeshCircle(elType,NL11,numofsegy,xmd);
    
    
end



if inclFilled
    [NL2,EL2] = squareRhombusInnerMesh(elType,numofsegy,x1,x2,x3,x4);
end

if (inclType == 'S' || inclType == 'R') && inclFilled
    %concatanering EL and NL
    NL = [NL1;NL2(4*numofsegy+1:end,:)];
    EL2 = EL2 + size(NL1,1)-4*numofsegy;
    EL = [EL1;EL2];

elseif (inclType == 'S' || inclType == 'R') && ~inclFilled
    NL = NL1;
    EL = EL1;
    
elseif inclType == 'C' && inclFilled
    NL = [NL1;NL11(4*numofsegy+1:end,:);NL2(4*numofsegy+1:end,:)];
    EL11 = EL11 + size(NL1,1)-4*numofsegy;
    EL2 = EL2 + size([NL1;NL11(4*numofsegy+1:end,:)],1)-4*numofsegy;
    EL = [EL1;EL11;EL2];
    
elseif inclType == 'C' && ~inclFilled
    NL = NL1;
    EL = EL1;
    

end




end

%############################################################
%############################################################
%############################################################

function EL = middleMesh(elType,NL,numofsegy,xmd)

numofsegc = 4*numofsegy;
pts = size(NL,1);
numNXpts = numofsegc;

if all(elType == 'D2QU4N')
    EL = zeros(pts,4);
    
    for i = 1:numofsegc %number of columns from the edges to the circle
        
        %for normal cases:
        EL1 = [(xmd+1)*(i-1) + 1 , (xmd+1)*(i) + 1 , (xmd+1)*(i) + 2 , (xmd+1)*(i-1) + 2];
        dummyEL = zeros(xmd+1,4);
        dummyEL(1,:) = EL1;
        
        %elements in btw
        for j = 1:xmd
            dummyEL(j+1,:) = dummyEL(j,:) + 1;
        end
        
        %the last element which its 2 nodes are touching the arc,
        dummyEL(end,3:4) = [pts-numNXpts+(i+1) pts-numNXpts+(i)];
        
        
        %last column
        if i == numofsegc
            
            EL1 = [(xmd+1)*(i-1) + 1 , 1 , 2 , (xmd+1)*(i-1) + 2];
            dummyEL = zeros(xmd+1,4);
            dummyEL(1,:) = EL1;
            
            %elements in btw
            for j = 1:xmd
                dummyEL(j+1,:) = dummyEL(j,:) + 1;
            end
            %the last element which its 2 nodes are touching the arc,
            dummyEL(end,3:4) = [pts-numNXpts+1, pts];
        end
        
        EL(1+(i-1)*(xmd+1):1+xmd+(i-1)*(xmd+1),:) = dummyEL;
        
    end
    
    
elseif all(elType == 'D2TR3N')
    
    EL = zeros(pts,3);
    for i = 1:numofsegc %number of columns from the edges to the circle
        
        %for normal cases:
        EL1 = [(xmd+1)*(i-1) + 1 , (xmd+1)*(i) + 1 , (xmd+1)*(i) + 2 , (xmd+1)*(i-1) + 2];
        dummyEL = zeros(xmd+1,4);
        dummyEL(1,:) = EL1;
        
        %elements in btw
        for j = 1:xmd
            dummyEL(j+1,:) = dummyEL(j,:) + 1;
        end
        
        %the last element which its 2 nodes are touching the arc,
        dummyEL(end,3:4) = [pts-numNXpts+(i+1) pts-numNXpts+(i)];
        
        
        %last column
        if i == numofsegc
            
            EL1 = [(xmd+1)*(i-1) + 1 , 1 , 2 , (xmd+1)*(i-1) + 2];
            dummyEL = zeros(xmd+1,4);
            dummyEL(1,:) = EL1;
            
            %elements in btw
            for j = 1:xmd
                dummyEL(j+1,:) = dummyEL(j,:) + 1;
            end
            %the last element which its 2 nodes are touching the arc,
            dummyEL(end,3:4) = [pts-numNXpts+1, pts];
        end
        
        dummyELtri = zeros(2*size(dummyEL,1),3);
        
        if ( (i>=1)&&(i<=numofsegy) ) || ( (i>=1 + 2*numofsegy)&&(i<=1 + 3*numofsegy) )
            dummyELtri(1:size(dummyEL,1),:) = dummyEL(:,[1 3 4]);
            dummyELtri(1+size(dummyEL,1):end,:) = dummyEL(:,[1 2 3]);
        else
            
            dummyELtri(1:size(dummyEL,1),:) = dummyEL(:,[1 2 4]);
            dummyELtri(1+size(dummyEL,1):end,:) = dummyEL(:,[2 3 4]);
        end
        
        EL(1+(i-1)*2*(xmd+1):(i)*2*(xmd+1),:) = dummyELtri;
        
    end
    
end
%removing zero rows from the EL
EL(EL(:,1)==0,:) = [] ;





end




%############################################################
%############################################################
%############################################################

function EL = middleMeshCircle(elType,NL,numofsegy,xmd)

numofsegc = 4*numofsegy;
pts = size(NL,1);
numNXpts = numofsegc;

if all(elType == 'D2QU4N')
    EL = zeros(pts,4);
    
    for i = 1:numofsegc %number of columns from the edges to the circle
        
        %for normal cases:
        EL1 = [(i-1) + 1 , (i-1) + 2,numofsegc+1+(i)*(xmd),numofsegc+1+(i-1)*(xmd)];
        dummyEL = zeros(xmd+1,4);
        dummyEL(1,:) = EL1;
        
        %elements in btw
        for j = 1:xmd
            dummyEL(j+1,:)=[dummyEL(j,4),dummyEL(j,3),numofsegc+j+1+(i)*(xmd),numofsegc+j+1+(i-1)*(xmd)];
        end
        
        %the last element which its 2 nodes are touching the arc,
        dummyEL(end,3:4) = [pts-numNXpts+(i+1) pts-numNXpts+(i)];
        
        
        %last column
        if i == numofsegc
            
            EL1 = [(i-1) + 1 , 1 , numofsegc+1  , numofsegc+1+(i-1)*(xmd)];
            dummyEL = zeros(xmd+1,4);
            dummyEL(1,:) = EL1;
            
            %elements in btw
            for j = 1:xmd
                dummyEL(j+1,:)=[dummyEL(j,4),dummyEL(j,3),numofsegc+1+j,numofsegc+j+1+(i-1)*(xmd)];
            end
            %the last element which its 2 nodes are touching the arc,
            dummyEL(end,3:4) = [pts-numNXpts+1, pts];
        end
        
        EL(1+(i-1)*(xmd+1):1+xmd+(i-1)*(xmd+1),:) = dummyEL;
        
    end
    
    
elseif all(elType == 'D2TR3N')
    
    EL = zeros(pts,3);
    for i = 1:numofsegc %number of columns from the edges to the circle
        
        %for normal cases:
        EL1 = [(i-1) + 1 , (i-1) + 2,numofsegc+1+(i)*(xmd),numofsegc+1+(i-1)*(xmd)];
        dummyEL = zeros(xmd+1,4);
        dummyEL(1,:) = EL1;
        
        %elements in btw
        for j = 1:xmd
            dummyEL(j+1,:)=[dummyEL(j,4),dummyEL(j,3),numofsegc+j+1+(i)*(xmd),numofsegc+j+1+(i-1)*(xmd)];
        end
        
        %the last element which its 2 nodes are touching the arc,
        dummyEL(end,3:4) = [pts-numNXpts+(i+1) pts-numNXpts+(i)];
        
        
        %last column
        if i == numofsegc
            
            EL1 = [(i-1) + 1 , 1 , numofsegc+1  , numofsegc+1+(i-1)*(xmd)];
            dummyEL = zeros(xmd+1,4);
            dummyEL(1,:) = EL1;
            
            %elements in btw
            for j = 1:xmd
                dummyEL(j+1,:)=[dummyEL(j,4),dummyEL(j,3),numofsegc+1+j,numofsegc+j+1+(i-1)*(xmd)];
            end
            %the last element which its 2 nodes are touching the arc,
            dummyEL(end,3:4) = [pts-numNXpts+1, pts];
        end
        
        dummyELtri = zeros(2*size(dummyEL,1),3);
        
        if ( (i>=1)&&(i<=numofsegy) ) || ( (i>=1 + 2*numofsegy)&&(i<=1 + 3*numofsegy) )
            dummyELtri(1:size(dummyEL,1),:) = dummyEL(:,[1 3 4]);
            dummyELtri(1+size(dummyEL,1):end,:) = dummyEL(:,[1 2 3]);
        else
            
            dummyELtri(1:size(dummyEL,1),:) = dummyEL(:,[1 2 4]);
            dummyELtri(1+size(dummyEL,1):end,:) = dummyEL(:,[2 3 4]);
        end
        
        EL(1+(i-1)*2*(xmd+1):(i)*2*(xmd+1),:) = dummyELtri;
        
    end
    
end
%removing zero rows from the EL
EL(EL(:,1)==0,:) = [] ;





end






%############################################################
%############################################################
%############################################################



function [NLnew,ELnew] = squareRhombusInnerMesh(elType,numofsegy,x1,x2,x3,x4)

NLouter = zeros(4*numofsegy,2);

dummy1 = lineGen(x1,x2,numofsegy);
NLouter(1:numofsegy,:) = dummy1(1:end-1,:);
% dummy1(1,:)=[];
% dummy1(end,:)=[];

dummy2 = lineGen(x2,x3,numofsegy);
NLouter(numofsegy+1:2*numofsegy,:) = dummy2(1:end-1,:);
dummy2(1,:)=[];
dummy2(end,:)=[];

dummy3 = lineGen(x3,x4,numofsegy);
NLouter(2*numofsegy+1:3*numofsegy,:) = dummy3(1:end-1,:);
% dummy3(1,:)=[];
% dummy3(end,:)=[];

dummy4 = lineGen(x4,x1,numofsegy);
NLouter(3*numofsegy+1:4*numofsegy,:) = dummy4(1:end-1,:);
dummy4(1,:)=[];
dummy4(end,:)=[];
dummy4(1:1:end,:) = dummy4(end:-1:1,:);

% x11 = [dummy1(1,1) dummy4(end,2)];
% x22 = [dummy1(end,1) x11(2)];
% x33 = [x22(1) dummy2(end,2)];
% x44 = [x11(1) x33(2)];

% lowerLine = lineGen(x11,x22,numofsegy-2);
% upperLine = lineGen(x11,x44,numofsegy-2);

segNum = numofsegy-2;

NL1 = zeros((segNum+1)^2,2); %innerCore
cnt = 1;

for i = 1:size(dummy4,1)
%     dummy = lineGen(upperLine(i,:),[x22(1) upperLine(i,2)],numofsegy-2);
    dummy = lineGen(dummy4(i,:),dummy2(i,:),numofsegy);
    dummy(1,:) = [];
    dummy(end,:) = [];
    
    NL1(cnt:cnt+size(dummy,1)-1,:) = dummy;
    cnt = cnt + size(dummy,1);
    
end




if all(elType == 'D2QU4N')
    
    EL1 = zeros(segNum^2,4);

    cnt = 1;
    %CRETING THE INNER CORE
    for i = 1:segNum
        dummyEL = [1+(i-1)*(segNum+1),2+(i-1)*(segNum+1),2+(i)*(segNum+1),1+(i)*(segNum+1)];

        v1 = dummyEL(1):dummyEL(1) + segNum - 1;
        v2 = dummyEL(2):dummyEL(2) + segNum - 1;
        v3 = dummyEL(3):dummyEL(3) + segNum - 1;
        v4 = dummyEL(4):dummyEL(4) + segNum - 1;
        dummyEL1 = [v1' v2' v3' v4'];
        
        EL1(cnt:cnt+size(dummyEL1,1)-1,:) = dummyEL1;
        cnt = cnt + size(dummyEL1,1);
    end
    
    %CREATING THE OUTER SHELL
    EL0 = zeros(4*(numofsegy-1),4);
    for i = 1:4
        
        if i == 1
            EL0(1,:) = [1 2 4*numofsegy+1 4*numofsegy];
            col1 = 2:numofsegy-1;
            col2 = 3:numofsegy;
            col3 = 4*numofsegy+2:5*numofsegy-1;
            col4 = 4*numofsegy+1:5*numofsegy-2;
            EL0(2:numofsegy-1,:) = [col1' col2' col3' col4'];

        elseif i == 2
            EL0(numofsegy,:) = [numofsegy,numofsegy+1,numofsegy+2,5*numofsegy-1];
            col1 = 5*numofsegy-1:numofsegy-1:(numofsegy-1)^2 + (4*numofsegy) - (numofsegy-1);
            col2 = numofsegy+2:2*numofsegy-1;
            col3 = col2+1;
            col4 = col1 + numofsegy-1;
            EL0(numofsegy+1:2*(numofsegy-1),:) = [col1' col2' col3' col4'];
            
        elseif i == 3
            EL0(2*(numofsegy-1)+1,:) = [col4(end) col3(end) 1+2*numofsegy 2+2*numofsegy];
            col1 = EL0(2*(numofsegy-1)+1,1)-1:-1:(4*numofsegy+1) +...
                        (numofsegy-1)*(numofsegy-2);
            col2 = col1 + 1;
            col3 = 2+2*numofsegy:3*numofsegy-1;
            col4 = col3+1;
            EL0(2*(numofsegy-1)+2:3*(numofsegy-1),:) = [col1' col2' col3' col4'];
            
        elseif i == 4
            EL0(3*(numofsegy-1)+1,:) = [3*numofsegy+2 col1(end) 3*numofsegy 3*numofsegy+1];
            col1 = 3*(numofsegy+1):4*numofsegy;
            col2 = EL0(3*(numofsegy-1)+1,2)-(numofsegy-1):1-numofsegy:4*numofsegy+1;
            col3 = col2 + numofsegy - 1;
            col4 = col1 - 1;
            EL0(3*(numofsegy-1)+2:end,:) = [col1' col2' col3' col4'];
        end

    end
    
    %putting it all together
    EL1= EL1 + 4*numofsegy;
    
    NLnew = [NLouter;NL1];
    ELnew = [EL0;EL1];
  
    
elseif all(elType == 'D2TR3N')

    EL1 = zeros(2*segNum^2,3);

    cnt = 1;
    %CRETING THE INNER CORE
    for i = 1:segNum
        dummyEL1 = [1+(i-1)*(segNum+1),2+(i-1)*(segNum+1),2+(i)*(segNum+1),1+(i)*(segNum+1)];
        dummyEL2 = [1+(i-1)*(segNum+1),2+(i-1)*(segNum+1),2+(i)*(segNum+1),1+(i)*(segNum+1)];

        v1 = dummyEL1(1):dummyEL1(1) + segNum - 1;
        v2 = dummyEL1(2):dummyEL1(2) + segNum - 1;
        v3 = dummyEL1(3):dummyEL1(3) + segNum - 1;
        dummyEL1 = [v1' v2' v3'];

        v1 = dummyEL2(1):dummyEL2(1) + segNum - 1;
        v3 = dummyEL2(3):dummyEL2(3) + segNum - 1;
        v4 = dummyEL2(4):dummyEL2(4) + segNum - 1;
        dummyEL2 = [v1' v3' v4'];

        
        EL1(cnt:cnt+size(dummyEL1,1)-1,:) = dummyEL1;
        cnt = cnt + size(dummyEL1,1);

        EL1(cnt:cnt+size(dummyEL2,1)-1,:) = dummyEL2;
        cnt = cnt + size(dummyEL2,1);
    end

    cnt = 1;
    %CREATING THE OUTER SHELL
    EL0 = zeros(8*(numofsegy-1),3);
    for i = 1:4
        
        if i == 1
            EL0(1,:) = [1 2 4*numofsegy+1];
            EL0(2,:) = [1 4*numofsegy+1 4*numofsegy];
            cnt = cnt + 2;

            col1 = 2:numofsegy-1;
            col2 = 3:numofsegy;
            col3 = 4*numofsegy+2:5*numofsegy-1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col2' col3'];
            cnt = cnt+length(col1);

            col1 = 2:numofsegy-1;
            col3 = 4*numofsegy+2:5*numofsegy-1;
            col4 = 4*numofsegy+1:5*numofsegy-2;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col3' col4'];
            cnt = cnt+length(col1);

        elseif i == 2
            EL0(cnt,:) = [numofsegy,numofsegy+1,numofsegy+2];
            cnt = cnt + 1;
            EL0(cnt,:) = [numofsegy,numofsegy+2,5*numofsegy-1];
            cnt = cnt + 1;

            col1 = 5*numofsegy-1:numofsegy-1:(numofsegy-1)^2 + (4*numofsegy) - (numofsegy-1);
            col2 = numofsegy+2:2*numofsegy-1;
            col3 = col2+1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col2' col3'];
            cnt = cnt+length(col1);

            col1 = 5*numofsegy-1:numofsegy-1:(numofsegy-1)^2 + (4*numofsegy) - (numofsegy-1);
            col3 = col2+1;
            col4 = col1 + numofsegy-1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col3' col4'];
            cnt = cnt+length(col1);
            
        elseif i == 3
            EL0(cnt,:) = [col4(end) col3(end) 1+2*numofsegy];   
            cnt = cnt + 1;
            EL0(cnt,:) = [col4(end) 1+2*numofsegy 2+2*numofsegy];
            cnt = cnt + 1;

            col1 = 4*numofsegy+((numofsegy-1)^2)-1:-1:(4*numofsegy+1)+(numofsegy-1)*(numofsegy-2);
            col2 = col1 + 1;
            col3 = 2+2*numofsegy:3*numofsegy-1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col2' col3'];
            cnt = cnt+length(col1);

            col1 = 4*numofsegy+((numofsegy-1)^2)-1:-1:(4*numofsegy+1)+(numofsegy-1)*(numofsegy-2);
            col3 = 2+2*numofsegy:3*numofsegy-1;
            col4 = col3+1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col3' col4'];
            cnt = cnt+length(col1);
            
        elseif i == 4
            EL0(cnt,:) = [3*numofsegy+2 col1(end) 3*numofsegy];
            cnt = cnt + 1;
            EL0(cnt,:) = [3*numofsegy+2 3*numofsegy 3*numofsegy+1];
            cnt = cnt + 1;

            col1 = 3*(numofsegy+1):4*numofsegy;
            col2 = 4*numofsegy+((numofsegy-1)^2)-(numofsegy-2)-(numofsegy-1):1-numofsegy:4*numofsegy+1;
            col3 = col2 + numofsegy - 1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col2' col3'];
            cnt = cnt+length(col1);

            col1 = 3*(numofsegy+1):4*numofsegy;
            col3 = col2 + numofsegy - 1;
            col4 = col1 - 1;
            EL0(cnt:cnt+length(col1)-1,:) = [col1' col3' col4'];
            cnt = cnt+length(col1);
        end

    end
    
    %putting it all together
    EL1= EL1 + 4*numofsegy;
    
    NLnew = [NLouter;NL1];
    ELnew = [EL0;EL1];
    
end


end




%############################################################
%############################################################
%############################################################



function linePts = lineGen(x1,x2,n)
%Generates a line with n segments btw coordinates x1 and x2

linePts = [( linspace(x1(1),x2(1),n+1) )' ( linspace(x1(2),x2(2),n+1) )'];

end
