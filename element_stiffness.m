function k = element_stiffness(nl,NL)
x = NL(nl,:);
NPE = size(x,1);
PD = size(x,2);
GPE = 1;

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
                            
                            S(a,c) = S(a,c) + grad(b,i)*constitutive(a,b,c,d)*grad(d,j)*det(J)*alpha;
                            
                            
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


end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [xi,eta,alpha] = GaussPoints(NPE,GPE,gp)

if (NPE == 3)
    
    xi = 1/3;
    eta = 1/3;
    alpha = 1*0.5; %triangular
    
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
        
        
    end
    
end



end

function C = constitutive (i,j,k,l)

%i = a
%j = b
%k = c
%l = d

E = 80/3;
nu = 1/3;
C = (E/(2*(1+nu)))*((i==l)*(j==k) + (i==k)*(j==l))...
    + ((E*nu)/(1-nu^2))*(i==j)*(k==l);

end