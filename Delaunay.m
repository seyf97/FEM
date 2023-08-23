function EL = Delaunay(NL,numNXpts,edgepts)
%                           DeLaunay Triangulation algorithm




%edgepts are the edge coordinates of the whole shape.
NoN = size(NL,1);
maxelnum = 4*NoN;
EL = zeros(maxelnum,3);

%Void Nodes that will be avoided
NX = zeros(numNXpts,1);
%Writing from left to right, the top boundary
NX(1:numNXpts,1) = size(NL,1):-1:size(NL,1)-(numNXpts)+1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%counter for element list
cnt = 1;

%Number of neighbouring pts to be searched for a triangle (12 default)
evalpts = 12;

if size(NL,1) <= evalpts
    evalpts = size(NL,1) - 1;
end

for i = 1:NoN
    
    
    
    disti = sqrt((NL(i,1)-NL(:,1)).^2 + (NL(i,2)-NL(:,2)).^2);
    
    [distances,index] = sort(disti);
    nl = index(2:evalpts+1);
    
    for j = 1:evalpts-1
        
        for k = 1:evalpts
            
            %1. Making sure not to count the same points in the 12 pt
            %cluster
            
            if j>=k
                continue
            end
            
            %2.Checking if the 3 points have already formed and element
            if any(sum( (( EL == i )+( EL == nl(j) )+( EL == nl(k) )),2) == 3 )
                continue
            end
            
% %                             Test pts reserved for debugging

%             testtt  = [50 1 24];
%             if sum(any(i == testtt) + any(nl(j) == testtt)...
%                     + any(nl(k) == testtt)) == 3
%                 
%                 i;
%                 
%             end
            
            
            %3.checking whether ALL 3 of the pointers are on the forbidden nodes, NX
            if any(sum( (( NX == i )+( NX == nl(j) )+( NX == nl(k) )),1) == 3 )...
                    && ~rightRay(edgepts,sum(NL([i nl(j) nl(k)],:),1)./3)
                continue
            end
            
            %4.Checking if the 3 pts are a line
            x = [NL(i,:);NL(nl(j),:);NL(nl(k),:)];
            v1 = x(2,:) - x(1,:);
            v2 = x(3,:) - x(2,:);
            v = [v1;v2];
            %if the cross product goes to zero, then skip the 3 lines. The
            %reason for putting the value in between +-1.e-10 is because
            %the value converges to 0, but matlab recognizes is as
            %infitesimally small due to decimal points
            if (v(1,1)*v(2,2) - v(2,1)*v(1,2))<1.e-10 && (v(1,1)*v(2,2) - v(2,1)*v(1,2))>-1.e-10
                continue
            end
            
            
            %Creating a circle with the 3 points i j k
            [xc,yc,r] = circle(x);
            %rounding in order to prevent numerical errors
            r = round(r,7);
            
            XXdummy = NL;
            %deleting the three points we are on
            XXdummy([i nl(j) nl(k)],:) = [];
            
            %calculating the distance from the center of the circumcircle
            %to the other points
            dist = round(sqrt((xc-XXdummy(:,1)).^2 + (yc-XXdummy(:,2)).^2),7);
            
            %checking if all the other points are not inside the circle
            %if yes, then arranging the nodes so that they are
            %counter-clockwise and putting them in the EL
            if all(dist >= r)
                
                %Checking if the triangle formed out of the 3 pts pass through
                %NX boundaries i nl(j) nl(k)
                %Checking for all Boundaries
                
                [list,indices] = polygon(NL,xc,yc,r,EL,edgepts);
                
                %putting the elements into the elements list
                EL(cnt:cnt+size(list,1)-1,:) = list;
                
                %putting the nodes to the list to avoid, since we don't
                %want to overlap triangles
                cnt = cnt + size(list,1);
            end
            
            
            
        end
        
    end
    
end


%removing zero rows from the EL
EL(EL(:,1)==0,:) = [] ;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xc,yc,r] = circle(x)

a1 = (x(1,1) - x(2,1))/(x(2,2) - x(1,2));
a2 = (x(2,1) - x(3,1))/(x(3,2) - x(2,2));

ym1 = (x(1,2) + x(2,2))/2;
ym2 = (x(2,2) + x(3,2))/2;

xm1 = (x(1,1) + x(2,1))/2;
xm2 = (x(2,1) + x(3,1))/2;

xc = (ym2 - ym1 + a1*xm1 - a2*xm2)/(a1-a2);
yc = ym1 + a1*xc - a1*xm1;



if round((x(2,2) - x(1,2)),7) == 0
    xc = xm1;
    yc = ym2 + a2*xc - a2*xm2;
end

if round((x(3,2) - x(2,2)),7) == 0
    xc = xm2;
    yc = ym1 + a1*xc - a1*xm1;
end


A = sqrt( (x(1,1) - x(2,1))^2 + (x(1,2)-x(2,2))^2 );
B = sqrt( (x(2,1) - x(3,1))^2 + (x(2,2)-x(3,2))^2 );
C = sqrt( (x(3,1) - x(1,1))^2 + (x(3,2)-x(1,2))^2 );

alpha = acos( (B^2 + C^2 - A^2)/(2*B*C) );
r = A/(2*sin(alpha));


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [list,indices] = polygon(NL,xc,yc,r,EL,edgepts)

edgeFlag = false;
flagEl = false;
%calculating the distance from the center of the circumcircle
%to the other points

%distance of all the points to the current circle
distPoly = sqrt((xc-NL(:,1)).^2 + (yc-NL(:,2)).^2);

%rounding to get rid of decimal errors
distPoly = round(distPoly,7);
r = round(r,7);

%counting the number of vertices of the shape
numofV = sum(r==distPoly);

%initializing the number of elements that will be produced
list = zeros(numofV-2,3);

%indices of the nodes that will be put inside
indices = zeros(1,numofV);

%corresponding angles btw the +x axis
angles = indices;

%finding the nodes
dummy1 = (distPoly == r);
dummy2 = 1:length(distPoly);
dummy3 = dummy1.*(dummy2)';
dummy3(dummy3==0)=[];
indices = dummy3;


for it = 1:length(indices)
    
    %finding the angle of the corresponding index
    s1 = round((NL(indices(it),2) - yc)/r,5);
    c1 = round((NL(indices(it),1) - xc)/r,5);
    if s1>0 && c1>0 %1st quadrant
        angle1 = atand(s1/c1);
    elseif s1>0 && c1<0 %2nd quadrant
        angle1 = 180 - atand(-s1/c1);
    elseif s1<0 && c1<0 %3rd quadrant
        angle1 = 180 + atand(s1/c1);
    elseif s1<0 && c1>0 %4th quadrant
        angle1 = 360 - atand(-s1/c1);
    elseif s1 == 0 && c1 == 1 % 0 deg
        angle1 = 0;
    elseif s1 == 1 && c1 == 0 %90 deg
        angle1 = 90;
    elseif s1 == 0 && c1 == -1 % 180 deg
        angle1 = 180;
    elseif s1 == -1 && c1 == 0 % 270 deg
        angle1 = 270;
    end
    
    angles(it) = angle1;
    
    %looking for a point which is on the shape's edge
    if any(sum(round(NL(indices(it),:),7) == round(edgepts,7),2) == 2)
        edgeFlag = true;        
    end
    
end


[sortedAng,indicesOfAng] = sort(angles);


for it = 1:size(list,1)
    
    %1st node on the triangle = indices(indicesofAng(1))
    %creating the triangles counter clockwise
    list(it,:) = [indices(indicesOfAng(1)),indices(indicesOfAng(1+it)),indices(indicesOfAng(2+it))];
    
    %Checking if the element is outside the boundaries
    if ~rightRay(edgepts,sum(NL([list(it,:)],:),1)./3)
        list(it,:) = [0 0 0];
        continue
    end
    
    %Checking if the 3 points have already formed and element
    if any(sum( (( EL == list(it,1) )+( EL == list(it,2) )+( EL == list(it,3) )),2) == 3 )
        list(it,:) = [0 0 0];
        flagEl = true;
        
    else
        num1 = list(it,1);
        num2 = list(it,2);
        num3 = list(it,3);
        
        %adding the last edgept to form the last line
        edgeptsNew = [edgepts;edgepts(1,:)];
        flagContinue = false;
        
        %Checking whether the element sides cross boundaries
        for pp = 1:size(edgepts,1)
            
            [flag1, x1] = lineIntersection(NL(num1,:),NL(num2,:),edgeptsNew(pp,:)...
                ,edgeptsNew(pp+1,:) );
            [flag2, x2] = lineIntersection(NL(num2,:),NL(num3,:),edgeptsNew(pp,:)...
                ,edgeptsNew(pp+1,:) );
            [flag3, x3] = lineIntersection(NL(num1,:),NL(num3,:),edgeptsNew(pp,:)...
                ,edgeptsNew(pp+1,:) );
            
            if flag1 == 2  || flag2 == 2 || flag3 == 2 %if any of the lines pass thru boundary
                flagContinue = true;
                break
            end
            
        end
        
        if flagContinue
            list(it,:) = [0 0 0];
        end
        
        
    end
    
    
end

%erasing the all zero rows which we didnt want to include
list(list(:,1)==0,:) = [] ;

%IF A PT IS ON THE EDGE AND THE ELEMENT WAS MISSED
if edgeFlag && ~flagEl && isempty(list)
    
    %initializing the number of elements that will be produced
    list = zeros(numofV-2,3);
    
    indicesOfAng = [indicesOfAng indicesOfAng(1)];
    indicesOfAng(1) = [];
    
    for it = 1:size(list,1)
        
        %1st node on the triangle = indices(indicesofAng(1))
        %creating the triangles counter clockwise
        list(it,:) = [indices(indicesOfAng(1)),indices(indicesOfAng(1+it)),indices(indicesOfAng(2+it))];
        
        %Checking if the element is outside the boundaries
        if ~rightRay(edgepts,sum(NL([list(it,:)],:),1)./3)
            list(it,:) = [0 0 0];
            continue
        end
        
        %Checking if the 3 points have already formed and element
        if any(sum( (( EL == list(it,1) )+( EL == list(it,2) )+( EL == list(it,3) )),2) == 3 )
            list(it,:) = [0 0 0];
            
        else
            num1 = list(it,1);
            num2 = list(it,2);
            num3 = list(it,3);
            
            %adding the last edgept to form the last line
            edgeptsNew = [edgepts;edgepts(1,:)];
            flagContinue = false;
            
            %Checking whether the element sides cross boundaries
            for pp = 1:size(edgepts,1)
                
                [flag1, x1] = lineIntersection(NL(num1,:),NL(num2,:),edgeptsNew(pp,:)...
                    ,edgeptsNew(pp+1,:) );
                [flag2, x2] = lineIntersection(NL(num2,:),NL(num3,:),edgeptsNew(pp,:)...
                    ,edgeptsNew(pp+1,:) );
                [flag3, x3] = lineIntersection(NL(num1,:),NL(num3,:),edgeptsNew(pp,:)...
                    ,edgeptsNew(pp+1,:) );
                
                if flag1 == 2  || flag2 == 2 || flag3 == 2 %if any of the lines pass thru boundary
                    flagContinue = true;
                    break
                end
                
            end
            
            if flagContinue
                list(it,:) = [0 0 0];
            end
            
            
        end
        
        
    end
    
    
    
    %erasing the all zero rows which we didnt want to include
    list(list(:,1)==0,:) = [] ;
end



end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function flagIn = rightRay(edgePts,x)
flagIn = false;
flagPt = false;

%takes in an x point and edge points and draws a ray to the right
%if the pt lies inside, odd number of intersections

%if the ray passes edge points 2*n times, or the lines 2*n times,
%it is outside the boundaries

%if a point is ON the edgePts or lie ON the lines, it will be considered
%outside

%if a point is on the edgePts, exit function

if any(sum(x == edgePts,2) == 2)
    flagIn = false;
    return
end

if any((x == edgePts(:,2)) == [0 1])
    flagPt = true;
end


%modifying so that we also obtain the last line in the for loop
newEdgePts = [edgePts;edgePts(1,:)];

if flagPt
    dummy2 = x;
    cnt1 = 0;
    cnt2 = 0;
    for k = 1:2
        
        %if the ray passes only thru edge pts, repeat process by shifting x upwards
        %0.0001 units and downwards 0.0001 units
        if k == 1
            x = dummy2 + [0 0.0001];
        else
            x = dummy2 + [0 -0.0001];
        end
        
        %for lines to the right of x
        cnt = 0;
        for i = 1:size(edgePts,1)
            
            %checks if the line is on the left of x
            if   x(1) > max(newEdgePts(i:i+1,1))
                continue
            end
            
            %checks if the line passes on a horizontal boundary
            if x(2) == newEdgePts(i,2) && x(2) == newEdgePts(i+1,2)
                flagIn = false;
                return
            end
            
            %checks if the line is on a vertical boundary
            if x(1) == newEdgePts(i,1) && x(1) == newEdgePts(i+1,1)
                flagIn = false;
                return
            end
            
            
            %creating the line equation (xline,yline)
            if newEdgePts(i,1) == newEdgePts(i+1,1) %vertical line
                xline = newEdgePts(i,1);
            else
                mx = (newEdgePts(i,2) - newEdgePts(i+1,2))/(newEdgePts(i,1)-newEdgePts(i+1,1));
                nx = newEdgePts(i,2) - mx*(newEdgePts(i,1));
            end
            
            
            %checks if the ray passes through the line to the right
            if  x(2)<max(newEdgePts(i:i+1,2)) && x(2)>min(newEdgePts(i:i+1,2))
                
                %checks if  the pt is on the line (just for lines with m~=0 && m~=Inf)
                if round(mx*x(1) + nx,5) == round(x(2),5) && (newEdgePts(i,1) ~= newEdgePts(i+1,1))...
                        && (newEdgePts(i,2) ~= newEdgePts(i+1,2))
                    flagIn = false;
                    return
                end
                
                if (x(2) - nx)/mx > x(1) && (newEdgePts(i,1) ~= newEdgePts(i+1,1))
                    cnt = cnt + 1;
                elseif (newEdgePts(i,1) == newEdgePts(i+1,1))
                    cnt = cnt+1;
                end
                
                
            end
            
            
            
        end
        
        if k == 1
            cnt1 = cnt;
            cnt = 0;
        else
            cnt2 = cnt;
            cnt = 0;
        end
        
        
    end
    
    if rem(cnt1,2) == 1 && rem(cnt2,2) == 1
        flagIn = true;
        return
    else
        flagIn = false;
        return
    end
    
else
    
    %for lines to the right of x
    cnt = 0;
    for i = 1:size(edgePts,1)
        
        %checks if the line is on the left of x
        if  x(1) > max(newEdgePts(i:i+1,1))
            continue
        end
        
        %checks if the line passes on a horizontal boundary
        if x(2) == newEdgePts(i,2) && x(2) == newEdgePts(i+1,2)
            flagIn = false;
            return
        end
        
        %checks if the line is on a vertical boundary
        if x(1) == newEdgePts(i,1) && x(1) == newEdgePts(i+1,1) && ...
                x(2)<max(newEdgePts(i:i+1,2)) && x(2)>min(newEdgePts(i:i+1,2))
            flagIn = false;
            return
        end
        
        
        %creating the line equation (xline,yline)
        if newEdgePts(i,1) == newEdgePts(i+1,1) %vertical line
            xline = newEdgePts(i,1);
        else
            mx = (newEdgePts(i,2) - newEdgePts(i+1,2))/(newEdgePts(i,1)-newEdgePts(i+1,1));
            nx = newEdgePts(i,2) - mx*(newEdgePts(i,1));
        end
        
        
        %checks if the ray passes through the line to the right
        if  x(2)<max(newEdgePts(i:i+1,2)) && x(2)>min(newEdgePts(i:i+1,2))
            
            %checks if  the pt is on the line (just for lines with m~=0 && m~=Inf)
            if round(mx*x(1) + nx,5) == round(x(2),5) && (newEdgePts(i,1) ~= newEdgePts(i+1,1))...
                    && (newEdgePts(i,2) ~= newEdgePts(i+1,2))
                flagIn = false;
                return
            end
            
            if (x(2) - nx)/mx > x(1) && (newEdgePts(i,1) ~= newEdgePts(i+1,1))
                cnt = cnt + 1;
            elseif (newEdgePts(i,1) == newEdgePts(i+1,1))
                cnt = cnt+1;
            end
            
        end
        
        
        
    end
    
    if rem(cnt,2) == 0
        flagIn = false;
    elseif rem(cnt,2) == 1
        flagIn = true;
    end
    
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [flag, x] = lineIntersection(q1,q2,z1,z2)
x = [0 0];

q1(1) = round(q1(1),7);
q1(2) = round(q1(2),7);

q2(1) = round(q2(1),7);
q2(2) = round(q2(2),7);

z1(1) = round(z1(1),7);
z1(2) = round(z1(2),7);

z2(1) = round(z2(1),7);
z2(2) = round(z2(2),7);

%Takes 2 lines and finds the intersection point if it exists
%First line is between q1 and q2, second line between z1 and z2

%If lines intersect inside the interval (x1,x2) --> x = [xcor ycord] flag=2
%If lines intersect at the boundaries x1||x2 --> x = [xcord ycord] flag=1

%If lines don't intersect --> x = [0 0], flag = -1;
%If lines Overlap --> x = [0 0], flag = -2;

%Finding the line  eq coefficients for both lines
% y = mx + n

%First checking if either of the lines are vertical

%rearranging the pts so that the first pt(q1) is on the left, and
%the second pt is on the right side of the space
if q1(1)>q2(1)
    dummy = q2;
    q2 = q1;
    q1 = dummy;
end

if z1(1)>z2(1)
    dummy = z2;
    z2 = z1;
    z1 = dummy;
end

if q1(1)>z1(1)
    x1 = q1;
else
    x1 = z1;
end

if q2(1)<z2(1)
    x2 = q2;
else
    x2 = z2;
end

%checking if the lower limit is larger than the upper limit,
%which results in the lines not sharing an interval
if x1(1)>x2(1)
    flag = -1;
    x = [0 0];
    return
end


%if q is a vertical line and z is not a vertical line
if q1(1) == q2(1) && z1(1) ~= z2(1)
    
    %finding line z,
    mz = (z2(2) - z1(2))/(z2(1) - z1(1));
    nz = z1(2) - mz*z1(1);
    y = mz*q1(1) + nz;
    y = round(y,7);
    
    %if y is btw the vertical line's y coordinates
    if  ( y<q2(2) && y>q1(2) ) || ( y<q1(2) && y>q2(2) )
        
        if z1(1)<q1(1) && q1(1)<z2(1)
            flag = 2;
            x = [q1(1) y];
            return
            
        elseif z1(1) == q1(1)
            flag = 1;
            x = [q1(1) y];
            return
            
        elseif z2(1) == q1(1)
            flag = 1;
            x = [q1(1) y];
            return
            
        end
        
        
        
    elseif y == q1(2)
        
        flag = 1;
        x = [q1(1) y];
        return
        
    elseif y == q2(2)
        
        flag = 1;
        x = [q1(1) y];
        return
        
    else
        flag = -1;
        x = [0 0];
        return
    end
    
    
    %if z is a vertical line and q is not a vertical line
elseif z1(1) == z2(1) && q1(1) ~= q2(1)
    
    %for line q: yq = mq(x) + nq
    mq = (q2(2) - q1(2))/(q2(1) - q1(1));
    nq = q1(2) - mq*q1(1);
    y = mq*z1(1) + nq;
    y = round(y,7);
    
    
    if  ( y<z2(2) && y>z1(2) ) || ( y<z1(2) && y>z2(2) )
        
        if q1(1)<z1(1) && z1(1)<q2(1)
            flag = 2;
            x = [z1(1) y];
            return
            
        elseif q1(1) == z1(1)
            flag = 1;
            x = [z1(1) y];
            return
            
        elseif q2(1) == z1(1)
            flag = 1;
            x = [z1(1) y];
            return
        end
        
        
    elseif y == z1(2)
        
        flag = 1;
        x = [z1(1) y];
        return
        
    elseif y == z2(2)
        
        flag = 1;
        x = [z1(1) y];
        return
    else
        flag = -1;
        x = [0 0];
        return
    end
    
    
    %both are a vertical line
elseif z1(1)== z2(1) && q1(1) == q2(1)
    
    if z1(1) ~= q1(1) %if x coordinates are not the same
        x = [0 0];
        flag = -1;
        return
    elseif z1(1) == q1(1)%if x coordinates are the same
        
        %are the lines in contact? if not, not intersecting
        if min([z1(2) z2(2)]) > max([q1(2) q2(2)]) ...
                || min([q1(2) q2(2)]) > max([z1(2) z2(2)])
            x = [0 0];
            flag = -1;
            return
            
            %are the max or min pts of the lines contacting?
            %this means they are on top of each other, and
            %touch each other at the middle
        elseif min([z1(2) z2(2)]) == max([q1(2) q2(2)]) ...
                || min([q1(2) q2(2)]) == max([z1(2) z2(2)])
            
            if min([z1(2) z2(2)]) == max([q1(2) q2(2)])
                y = min([z1(2) z2(2)]);
            else
                y = min([q1(2) q2(2)]);
            end
            
            x = [z1(1) y];
            flag = 1;
            return
            
            %the lines are in contact and they are not on top of each other
            %this leaves us with only one option: they are overlapping
        else
            x = [0 0];
            flag = -2;
        end
        
        
    end
    
else
    
    %for line q: yq = mq(x) + nq
    mq = (q2(2) - q1(2))/(q2(1) - q1(1));
    nq = q1(2) - mq*q1(1);
    
    %finding line z,
    mz = (z2(2) - z1(2))/(z2(1) - z1(1));
    nz = z1(2) - mz*z1(1);
    
    %yq - yz is the difference between the lines in the interval
    %we will check it for pts x1 and x2
    
    %if their difference does not change sign, no intersecting
    %if changes sign, intersects inside the interval
    %if the value is 0 at boundary, intersects at the interval bnd
    %if both 0, lines are overlapping
    %CHANGE THE BOTTOM AND UPPER VALUES SO THAT ITS IN THE INTERVAL WHEN TWO
    %LINES ARE DEFINED
    bottomValue = (mq-mz)*x1(1) + (nq-nz);
    upperValue = (mq-mz)*x2(1) + (nq-nz);
    
    bottomValue = round(bottomValue,5); %!!!ROUNDING IN ORDER TO ELIMINATE ERROR
    upperValue = round(upperValue,5);   %!!!ROUNDING IN ORDER TO ELIMINATE ERROR
    
    %if sign changes --> intersects in the interval
    if bottomValue*upperValue<0
        
        %finding the intersecting x coordinate
        x(1) = -(nq-nz)/(mq-mz);
        %putting in the x coordinate to either one of the line eq's
        x(2) = mq*x(1) + nq;
        %changing the flag
        flag = 2;
        return
        
    elseif bottomValue*upperValue>0 %no intersection
        flag = -1;
        x = [0 0];
        return
    elseif bottomValue == 0 && upperValue == 0 %overlap
        flag = -2;
        x = [0 0];
        return
    elseif bottomValue == 0 || upperValue == 0 %boundary intersect
        flag = 1;
        
        if bottomValue == 0
            x = [x1(1) x1(2)];
        else
            x = [x2(1) x2(2)];
        end
        return
        
    end
    
    
    
    
end



end