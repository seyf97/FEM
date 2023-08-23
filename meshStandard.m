function [EL,NL] = meshStandard(r,l,w,elType,numofsegy,xmd,cgx,cgy,inclType,inclFilled)

numofsegc = 4*(numofsegy);

if inclType == 'C'
    
    xc = r.*cosd(linspace(225,360+225,numofsegc+1));
    yc = r.*sind(linspace(225,360+225,numofsegc+1));
    
    
    %shifting the CG of the Circle
    xc = xc + cgx ;
    yc = yc + cgy;
    
    xc(end) = [];
    yc(end) = [];
    
    
    xc = [xc' yc'];
    
    
end




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

%adding points in between the wall and the void in the middle


%creating the Nodes List
NL = zeros(size(x0,1)*(1+xmd)+size(xc,1),2);

%The for loop below ensures that the boundary nodes are at the bottom of NL
for index = 1:size(x0,1)
    
    dummyx = linspace(x0(index,1),xc(index,1),xmd+2);
    dummyy = linspace(x0(index,2),xc(index,2),xmd+2);
    
    %since we want to put the boundary points to the bottom
    dummyx(end) = [];
    dummyy(end) = [];
    
    NL(1+(index-1)*(xmd+1):index*(xmd+1),:) = [dummyx' dummyy'];    
end
dummyIndex = size(NL,1) - size(xc,1) + 1;
NL(dummyIndex:end,:) = xc;


%putting the node coordinates in the NL
% NL = x;
%the boundary points are always at the bottom of the NL
numNXpts = numofsegc;



pts = size(NL,1);
maxelnum = pts*(pts-1)*(pts-2)/6;

%Nodes that will be avoided
%1st row: Boundaries
%n'th row: circle with more than 3 points on it
%!!! choosing 100 for the moment, the number will be updated.
NX = zeros(size(NL,1),100);
%Writing from left to right, the top boundary
NX(1,1:numNXpts) = size(NL,1):-1:size(NL,1)-(numNXpts)+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Creating quadrilateral elements
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
    
end



if all(elType == 'D2TR3N')
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

fprintf('Number of Elements: %d \nNumber of Nodes: %d \n',size(EL,1),size(NL,1));


%end of mesh function
end


