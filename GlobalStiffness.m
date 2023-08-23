function K = GlobalStiffness(EL,NL,ENL,GPE,Ei,nuI,Eo,nuO,numofsegy,xmd,NPE)

NoE = size(EL,1);
NoN = size(NL,1);
K = sparse(2*NoN,2*NoN);

for index = 1:NoE
    
    x = NL(EL(index,:),:);
    
    output = Stiffness(x,GPE,index,Ei,nuI,Eo,nuO,numofsegy,xmd);
    globalNodes = ENL(EL(index,:),7:8);
    dummy = zeros(1,2*NPE);
    dummy(1:2:end) = globalNodes(:,1);
    dummy(2:2:end) = globalNodes(:,2);
  
    K(dummy,dummy) = K(dummy,dummy) + output;
    
end




end



function output = Stiffness( x , GPE,index,Ei,nuI,Eo,nuO,numofsegy,xmd )

NPE = size(x,1);
PD = size(x,2);

coor = x';
k = zeros(PD*NPE,PD*NPE); %STIFFNESS IS X2, BECAUSE IN THERMAL ONLY SCALAR, BUT HERE X AND Y COMP.


for i = 1:NPE
    
    for j = 1:NPE
        
        K = zeros(PD,PD);
        S = zeros(PD,PD);
        
        for gp = 1:GPE
            
            J = zeros(PD,PD);
            
            grad = zeros(PD,NPE);
            
            [xi,eta,alpha] = GaussPoints(NPE,GPE,gp);
            
            grad_nat = grad_N_nat( NPE, xi, eta);
            
            J = coor*grad_nat';
            
            grad = inv(J)'*grad_nat;
            
            for a = 1:PD
                for c = 1:PD
                    for b = 1:PD
                        for d = 1:PD
                            %Now K is S
                            
                            S(a,c) = S(a,c) + grad(b,i)*constitutive(a,b,c,d,index,Ei,nuI,Eo,nuO,numofsegy,xmd,NPE)*grad(d,j)*det(J)*alpha;
                            
                            
                        end
                    end
                end
            end
                       
            
            
        end
        k((i-1)*PD+1:i*PD,(j-1)*PD+1:j*PD) = S;
        
    end
    
end

output = k;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = grad_N_nat(NPE, xi, eta) % gives with gaus pts implemented'

PD = 2;
result = zeros(PD,NPE);

if (NPE == 3)
    result(1,1) = 1;
    result(1,2) = 0;
    result(1,3) = -1;
    
    result(2,1) = 0;
    result(2,2) = 1;
    result(2,3) = -1;
    
elseif (NPE == 4)
    result(1,1) = -0.25*(1-eta);
    result(1,2) = +0.25*(1-eta);
    result(1,3) = +0.25*(1+eta);
    result(1,4) = -0.25*(1+eta);

    result(2,1) = -0.25*(1-xi);
    result(2,2) = -0.25*(1+xi);
    result(2,3) = 0.25*(1+xi);
    result(2,4) = 0.25*(1-xi);

   
elseif (NPE == 6)
    result(1,1) = -1 + 4*xi;
    result(1,2) = 0;
    result(1,3) = -3 + 4*(xi + eta);
    result(1,4) = 4*eta;
    result(1,5) = -4*eta;
    result(1,6) = -4*(-1+eta+2*xi);

    result(2,1) = 0;
    result(2,2) = -1 + 4*eta;
    result(2,3) = -3 + 4*(xi + eta);
    result(2,4) = 4*xi;
    result(2,5) = -4*(-1 + 2*eta + xi);
    result(2,6) = -4*xi;

   
elseif (NPE == 8)
    result(1,1) = +0.25*(1-eta)*(2*xi + eta);
    result(1,2) = +0.25*(1-eta)*(2*xi - eta);
    result(1,3) = +0.25*(1+eta)*(2*xi + eta);
    result(1,4) = +0.25*(1+eta)*(2*xi - eta);
    result(1,5) = -xi*(1-eta);
    result(1,6) = +0.5*(1-eta)*(1+eta);
    result(1,7) = -xi*(1+eta);
    result(1,8) = -0.5*(1-eta)*(1+eta);

    result(2,1) = +0.25*(1-xi)*(xi + 2*eta);
    result(2,2) = +0.25*(1+xi)*(-xi + 2*eta);
    result(2,3) = +0.25*(1+xi)*(xi + 2*eta);
    result(2,4) = +0.25*(1-xi)*(-xi + 2*eta);
    result(2,5) = -0.5*(1-xi)*(1+xi);
    result(2,6) = -(1+xi)*(eta);
    result(2,7) = 0.5*(1-xi)*(1+xi);
    result(2,8) = -(1-xi)*eta;


   
elseif (NPE == 9)
    result(1,1) = +0.25*(1-2*xi)*(1-eta)*eta;
    result(1,2) = -0.25*(1+2*xi)*(1-eta)*eta;
    result(1,3) = +0.25*(1+2*xi)*(1+eta)*eta;
    result(1,4) = -0.25*(1-2*xi)*(1+eta)*eta;
    result(1,5) = xi*eta*(1-eta);
    result(1,6) = +0.5*(1-eta)*(1+eta)*(1+2*xi);
    result(1,7) = -xi*eta*(1+eta);
    result(1,8) = -0.5*(1-eta)*(1+eta)*(1-2*xi);
    result(1,9) = -2*xi*(1-eta)*(1+eta);

    result(2,1) = +0.25*(1-xi)*xi*(1 - 2*eta);
    result(2,2) = -0.25*(1+xi)*xi*(1 - 2*eta);
    result(2,3) = +0.25*(1+xi)*xi*(1 + 2*eta);
    result(2,4) = -0.25*(1-xi)*xi*(1 + 2*eta);
    result(2,5) = 0.5*(1-xi)*(1+xi)*(2*eta - 1);
    result(2,6) = -(1+xi)*(eta)*xi;
    result(2,7) = 0.5*(1-xi)*(1+xi)*(1 + 2*eta);
    result(2,8) = (1-xi)*eta*xi;
    result(2,9) = -2*(1-xi)*(1+xi)*eta;


end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xi,eta,alpha] = GaussPoints(NPE,GPE,gp)

if (NPE == 3)
    
    if (GPE == 1)
        xi = 1/3;
        eta = 1/3;
        alpha = 1*0.5; %triangular
        
    elseif (GPE == 3)

        if gp == 1
            xi = 1/6;
            eta = 1/6;
            alpha = 1/3;
            
        elseif gp == 2
            xi = 4/6;
            eta = 1/6;
            alpha = 1/3;
            
        elseif gp == 3
            xi = 1/6;
            eta = 4/6;
            alpha = 1/3;
        end
        
    end
    
    
elseif (NPE == 6)
    
    if (GPE == 1)
        xi = 1/3;
        eta = 1/3;
        alpha = 1*0.5; %triangular
        
    elseif (GPE == 3)
        if gp == 1
            xi = 1/6;
            eta = 1/6;
            alpha = 1/3;
            
        elseif gp == 2
            xi = 4/6;
            eta = 1/6;
            alpha = 1/3;
            
        elseif gp == 3
            xi = 1/6;
            eta = 4/6;
            alpha = 1/3;
        end
        
    end
    
    
elseif (NPE == 4)
    
    if(GPE == 1)
        xi = 0;
        eta = 0;
        alpha = 4;
        
    elseif(GPE == 4)
        
        if gp == 1
            xi = -1/sqrt(3);
            eta = -1/sqrt(3);
            alpha = 1;
            
        elseif gp == 2
            xi = 1/sqrt(3);
            eta = -1/sqrt(3);
            alpha = 1;
            
        elseif gp == 3
            xi = 1/sqrt(3);
            eta = 1/sqrt(3);
            alpha = 1;
            
        elseif gp == 4
            xi = -1/sqrt(3);
            eta = 1/sqrt(3);
            alpha = 1;
        end

    elseif(GPE==9)   

        if gp == 1
            xi = -sqrt(3/5);
            eta = -sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 2
            xi = +sqrt(3/5);
            eta = -sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 3
            xi = +sqrt(3/5);
            eta = +sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 4
            xi = -sqrt(3/5);
            eta = +sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 5
            xi = 0;
            eta = -sqrt(3/5);
            alpha = 40/81;
            
        elseif gp == 6
            xi = sqrt(3/5);
            eta = 0;
            alpha = 40/81;
            
        elseif gp == 7
            xi = 0;
            eta = sqrt(3/5);
            alpha = 40/81;
            
        elseif gp == 8
            xi = -sqrt(3/5);
            eta = 0;
            alpha = 40/81;
            
        elseif gp == 9
            xi = 0;
            eta = 0;
            alpha = 64/81;
            
        end
        
    end
    
    
elseif (NPE == 8)
    
    if(GPE == 1)
        xi = 0;
        eta = 0;
        alpha = 4;
        
    elseif(GPE == 4)
        
        if gp == 1
            xi = -1/sqrt(3);
            eta = -1/sqrt(3);
            alpha = 1;
            
        elseif gp == 2
            xi = 1/sqrt(3);
            eta = -1/sqrt(3);
            alpha = 1;
            
        elseif gp == 3
            xi = 1/sqrt(3);
            eta = 1/sqrt(3);
            alpha = 1;
            
        elseif gp == 4
            xi = -1/sqrt(3);
            eta = 1/sqrt(3);
            alpha = 1;
        end
        
    elseif(GPE == 9)
        
        if gp == 1
            xi = -sqrt(3/5);
            eta = -sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 2
            xi = +sqrt(3/5);
            eta = -sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 3
            xi = +sqrt(3/5);
            eta = +sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 4
            xi = -sqrt(3/5);
            eta = +sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 5
            xi = 0;
            eta = -sqrt(3/5);
            alpha = 40/81;
            
        elseif gp == 6
            xi = sqrt(3/5);
            eta = 0;
            alpha = 40/81;
            
        elseif gp == 7
            xi = 0;
            eta = sqrt(3/5);
            alpha = 40/81;
            
        elseif gp == 8
            xi = -sqrt(3/5);
            eta = 0;
            alpha = 40/81;
            
        elseif gp == 9
            xi = 0;
            eta = 0;
            alpha = 64/81;
            
        end
        
    end
    
    
    
elseif (NPE == 9)
    
    if(GPE == 1)
        xi = 0;
        eta = 0;
        alpha = 4;
        
    elseif(GPE == 4)
        
        if gp == 1
            xi = -1/sqrt(3);
            eta = -1/sqrt(3);
            alpha = 1;
            
        elseif gp == 2
            xi = 1/sqrt(3);
            eta = -1/sqrt(3);
            alpha = 1;
            
        elseif gp == 3
            xi = 1/sqrt(3);
            eta = 1/sqrt(3);
            alpha = 1;
            
        elseif gp == 4
            xi = -1/sqrt(3);
            eta = 1/sqrt(3);
            alpha = 1;
        end
        
    elseif(GPE == 9)
        
        if gp == 1
            xi = -sqrt(3/5);
            eta = -sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 2
            xi = +sqrt(3/5);
            eta = -sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 3
            xi = +sqrt(3/5);
            eta = +sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 4
            xi = -sqrt(3/5);
            eta = +sqrt(3/5);
            alpha = 25/81;
            
        elseif gp == 5
            xi = 0;
            eta = -sqrt(3/5);
            alpha = 40/81;
            
        elseif gp == 6
            xi = sqrt(3/5);
            eta = 0;
            alpha = 40/81;
            
        elseif gp == 7
            xi = 0;
            eta = sqrt(3/5);
            alpha = 40/81;
            
        elseif gp == 8
            xi = -sqrt(3/5);
            eta = 0;
            alpha = 40/81;
            
        elseif gp == 9
            xi = 0;
            eta = 0;
            alpha = 64/81;
            
        end
        
    end
    
    
end



end

function C = constitutive(i,j,k,l,index,Ei,nuI,Eo,nuO,numofsegy,xmd,NPE)

%i = a
%j = b
%k = c
%l = d

if NPE == 3 || NPE == 6
    num = 2*numofsegy*(xmd+1)*4;
else
    num = numofsegy*(xmd+1)*4;
end


if index<=num
    E = Eo;
    nu = nuO;
else
    E = Ei;
    nu = nuI;
end

C = (E/(2*(1+nu)))*((i==l)*(j==k) + (i==k)*(j==l))...
    + ((E*nu)/(1-nu^2))*(i==j)*(k==l);

end