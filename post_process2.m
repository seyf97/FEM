function [NoN] = post_process2(NL,EL,ENL,Ei,nuI,Eo,nuO,numofsegy,xmd,GPE)

PD = 2;
NPE = size(EL, 2);
NoN = size(ENL, 1);
NoE = size(EL, 1);
scale = 1; %magnifies displacement/deformation of structure

[disp, stress, strain] = element_pp(NL, EL, ENL,Ei,nuI,Eo,nuO,numofsegy,xmd,GPE);

switch NPE
    case 3
        for i = 1:NoE
            %nl = EL(i, 1, NPE);
            nl = EL(i, 1:NPE);
            
            for j = 1:NPE
                X_1(j, i) = ENL(nl(j), 1) + scale*ENL(nl(j), 4*PD+1);
                Y_1(j, i) = ENL(nl(j), 2) + scale*ENL(nl(j), 4*PD+2);
            end
            
            for j = 1:NPE
                val = stress(i, :, 1, 1);
                stress_xx(j,i) = val;
            end
            
            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val;
            end
            
            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val;
            end
            
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val;
            end
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
            end
            for j = 1:NPE
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
            
            
            for j = 1:NPE
                val = strain(i, :, 1, 1);
                strain_xx(j,i) = val;
            end
            
            for j = 1:NPE
                val = strain(i, :, 1, 2);
                strain_xy(j,i) = val;
            end
            
            for j = 1:NPE
                val = strain(i, :, 2, 1);
                strain_yx(j,i) = val;
            end
            
            for j = 1:NPE
                val = strain(i, :, 2, 2);
                strain_yy(j,i) = val;
            end
        end
    case 4
        for i = 1:NoE
            %nl = EL(i, 1, NPE);
            nl = EL(i, 1:NPE);
            for j = 1:NPE
                X_1(j, i) = ENL(nl(j), 1) + scale*ENL(nl(j), 4*PD+1);
                Y_1(j, i) = ENL(nl(j), 2) + scale*ENL(nl(j), 4*PD+2);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 1);
                stress_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val(1, j);
            end
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
            end
            for j = 1:NPE
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
            
            for j = 1:NPE
                val = strain(i, :, 1, 1);
                strain_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 1, 2);
                strain_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 1);
                strain_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 2);
                strain_yy(j,i) = val(1, j);
            end
        end



    case 6
        for i = 1:NoE
            %nl = EL(i, 1, NPE);
            nl = EL(i, 1:NPE);
            for j = 1:NPE
                X_1(j, i) = ENL(nl(j), 1) + scale*ENL(nl(j), 4*PD+1);
                Y_1(j, i) = ENL(nl(j), 2) + scale*ENL(nl(j), 4*PD+2);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 1);
                stress_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val(1, j);
            end
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
            end
            for j = 1:NPE
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
            
            for j = 1:NPE
                val = strain(i, :, 1, 1);
                strain_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 1, 2);
                strain_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 1);
                strain_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 2);
                strain_yy(j,i) = val(1, j);
            end
        end





    case 8
        for i = 1:NoE
            %nl = EL(i, 1, NPE);
            nl = EL(i, 1:NPE);
            for j = 1:NPE
                X_1(j, i) = ENL(nl(j), 1) + scale*ENL(nl(j), 4*PD+1);
                Y_1(j, i) = ENL(nl(j), 2) + scale*ENL(nl(j), 4*PD+2);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 1);
                stress_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val(1, j);
            end
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
            end
            for j = 1:NPE
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
            
            for j = 1:NPE
                val = strain(i, :, 1, 1);
                strain_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 1, 2);
                strain_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 1);
                strain_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 2);
                strain_yy(j,i) = val(1, j);
            end
        end



    case 9
        for i = 1:NoE
            %nl = EL(i, 1, NPE);
            nl = EL(i, 1:NPE);
            for j = 1:NPE
                X_1(j, i) = ENL(nl(j), 1) + scale*ENL(nl(j), 4*PD+1);
                Y_1(j, i) = ENL(nl(j), 2) + scale*ENL(nl(j), 4*PD+2);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 1);
                stress_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 1, 2);
                stress_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 1);
                stress_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = stress(i, :, 2, 2);
                stress_yy(j,i) = val(1, j);
            end
            for j = 1:NPE
                disp_y(j, i) = abs(ENL(nl(j), 4*PD + 2));
            end
            for j = 1:NPE
                disp_x(j, i) = abs(ENL(nl(j), 4*PD + 1));
            end
            
            for j = 1:NPE
                val = strain(i, :, 1, 1);
                strain_xx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 1, 2);
                strain_xy(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 1);
                strain_yx(j,i) = val(1, j);
            end
            for j = 1:NPE
                val = strain(i, :, 2, 2);
                strain_yy(j,i) = val(1, j);
            end
        end
        
        
        
        
end

stressXX = stress(:,:,1,1);
stressYY = stress(:,:,2,2);
stressXY = stress(:,:,1,2);

strainXX = strain(:,:,1,1);
strainYY = strain(:,:,2,2);
strainXY = strain(:,:,1,2);

[avgStressXX,avgStressXY,avgStressYY,avgStrainXX,avgStrainXY,avgStrainYY] = ...
    findAvg(stressXX,stressYY,stressXY,strainXX,strainYY,strainXY,NoE,NoN,EL);





%  'EdgeColor','none'
figure('Color',[1 1 1]);
set(gca, 'color', [1 1 1])
hold on

subplot(3,3,1)
% figure(1)
patch(X_1, Y_1, avgStressXX, 'FaceColor' , 'interp');
title('\sigma_x_x')
colormap jet;
colorbar;
axis equal

subplot(3,3,2)
% figure(2)
patch(X_1, Y_1, avgStressXY, 'FaceColor' , 'interp');
title('\sigma_x_y')
colormap jet;
colorbar;
axis equal

subplot(3,3,3)
% figure(4)
patch(X_1, Y_1, avgStressYY, 'FaceColor' , 'interp');
title('\sigma_y_y')
colormap jet;
colorbar;
axis equal


% subplot(5,3,4)
% % figure(1)
% patch(X_1, Y_1, avgStressXX, 'FaceColor' , 'interp','EdgeColor','interp');
% title('Avg. \sigma_x_x')
% colormap jet;
% colorbar;
% % axis equal
% 
% subplot(5,3,5)
% % figure(2)
% patch(X_1, Y_1, avgStressXY, 'FaceColor' , 'interp','EdgeColor','interp');
% title('Avg. \sigma_x_y')
% colormap jet;
% colorbar;
% % axis equal
% 
% subplot(5,3,6)
% % figure(4)
% patch(X_1, Y_1, avgStressYY, 'FaceColor' , 'interp','EdgeColor','interp');
% title('Avg. \sigma_y_y')
% colormap jet;
% colorbar;
% % axis equal

subplot(3,3,4)
% figure(7)
patch(X_1, Y_1, avgStrainXX, 'FaceColor' , 'interp');
title('\epsilon_x_x')
colormap jet;
colorbar;
axis equal

subplot(3,3,5)
% figure(8)
patch(X_1, Y_1, avgStrainXY, 'FaceColor' , 'interp');
title('\epsilon_x_y')
colormap jet;
colorbar;
axis equal


subplot(3,3,6)
patch(X_1, Y_1, avgStrainYY, 'FaceColor' , 'interp');
title('\epsilon_y_y')
colormap jet;
colorbar;
axis equal

% subplot(5,3,10)
% % figure(7)
% patch(X_1, Y_1, avgStrainXX, 'FaceColor' , 'interp','EdgeColor','interp');
% title('Avg. \epsilon_x_x')
% colormap jet;
% colorbar;
% % axis equal
% 
% subplot(5,3,11)
% % figure(8)
% patch(X_1, Y_1, avgStrainXY, 'FaceColor' , 'interp','EdgeColor','interp');
% title('Avg. \epsilon_x_y')
% colormap jet;
% colorbar;
% % axis equal
% 
% 
% subplot(5,3,12)
% patch(X_1, Y_1, avgStrainYY, 'FaceColor' , 'interp','EdgeColor','interp');
% title('Avg. \epsilon_y_y')
% colormap jet;
% colorbar;
% % axis equal






subplot(3,3,7)
% figure(5)
patch(X_1, Y_1, disp_x, 'FaceColor' , 'interp');
title('\Delta X')
colormap jet;
colorbar;
axis equal

subplot(3,3,8)
% figure(6)
patch(X_1, Y_1, disp_y, 'FaceColor' , 'interp');
title('\Delta Y')
colormap jet;
colorbar;
axis equal


text = sprintf('Expansion Test \n %d TR Elements',size(EL,1));
sgtitle(text);





end







function [disp, stress, strain] = element_pp(NL, EL, ENL,Ei,nuI,Eo,nuO,numofsegy,xmd,GPE)

NoN = size(ENL, 1);
NoE = size(EL, 1);
NPE = size(EL, 2);
PD = size(NL, 2);

switch NPE
    case 3
        GPE = 1;
    case 4
        GPE = 4;
    case 6
        GPE = 3;
    case 8
        GPE = 4;
    case 9
        GPE = 9;
end

%we write displacements on the corners, stress and strain on gauss points
disp = zeros(NoE, NPE, PD, 1);
%look at contents of stress matrix
%we choose the element
%the element has a set number of gauss points (e.g. 4)
%each gauss point has a matrix of PD * PD
stress = zeros(NoE, GPE, PD, PD);
strain = zeros(NoE, GPE, PD, PD);

for e = 1:NoE
    
    nl = EL(e, 1:NPE); %node numbers connected to each element
    coor = zeros(PD, NPE);
    
    for i = 1:NPE
        for j = 1:PD
            coor(j, i) = NL(nl(i), j);
        end
    end
    
    u = zeros(PD, NoE);
    for i = 1:NPE
        for j = 1:PD
            u(j,i) = ENL(nl(i), 4*PD + j);
        end
    end
    
    for i = 1:NPE
        for j = 1:PD
            disp(e,i,j,1) = ENL(nl(i), 4*PD + j);
        end
    end
    
    for gp = 1:GPE
        epsilon = zeros(PD, PD);
        for i = 1:NPE
            J = zeros(PD, PD);
            grad = zeros(PD, NPE);
        end
    end
    
    
    for gp = 1:GPE
        
        epsilon = zeros(PD, PD);
        
        for i = 1:NPE
            
            J = zeros(PD, PD);
            grad = zeros(PD, NPE);
            [xi, eta, alpha] = GaussPoints(NPE, GPE, gp);
            grad_nat = grad_N_nat(NPE,xi,eta);
            J = coor*grad_nat';
            grad = inv(J)' * grad_nat;
            
            epsilon = epsilon + 1/2 * (dyad(u(:, i), grad(:, i)) + dyad(grad(:, i), u(:, i)));
        end
        
        sigma = zeros(PD, PD);
        for a = 1:PD
            for b = 1:PD
                for c = 1:PD
                    for d = 1:PD
                        sigma(a, b) = sigma(a, b) + constitutive(a,b,c,d,e,Ei,nuI,Eo,nuO,numofsegy,xmd,NPE) * epsilon(c, d);
                    end
                end
            end
        end
        
        for a = 1:PD
            for b = 1:PD
                strain(e,gp,a,b) = epsilon(a,b);
                stress(e,gp,a,b) = sigma(a,b);
            end
        end
        
    end
    
end


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


function A = dyad(u, v)
PD = 2;
A = zeros(PD, PD);
for i = 1:PD
    for j = 1:PD
        A(i, j) = u(i) * v(j);
    end
end


end



function [avgStressXX,avgStressXY,avgStressYY,avgStrainXX,avgStrainXY,avgStrainYY] = ...
    findAvg(stressXX,stressYY,stressXY,strainXX,strainYY,strainXY,NoE,NoN,EL)

avgStressXX = zeros(NoE,size(EL,2));
avgStressXY = zeros(NoE,size(EL,2));
avgStressYY = zeros(NoE,size(EL,2));
avgStrainXX = zeros(NoE,size(EL,2));
avgStrainXY = zeros(NoE,size(EL,2));
avgStrainYY = zeros(NoE,size(EL,2));

for i = 1:NoN
    
    %StressXX Avg
    n = sum(sum(i == EL));
    
    aStressXX = sum(sum((i==EL).*(stressXX)))/n;
    aStressXY = sum(sum((i==EL).*(stressXY)))/n;
    aStressYY = sum(sum((i==EL).*(stressYY)))/n;
    
    aStrainXX = sum(sum((i==EL).*(strainXX)))/n;
    aStrainXY = sum(sum((i==EL).*(strainXY)))/n;
    aStrainYY = sum(sum((i==EL).*(strainYY)))/n;
    
    avgStressXX = avgStressXX + aStressXX*(i==EL);
    avgStressXY = avgStressXY + aStressXY*(i==EL);
    avgStressYY = avgStressYY + aStressYY*(i==EL);
    
    avgStrainXX = avgStrainXX + aStrainXX*(i==EL);
    avgStrainXY = avgStrainXY + aStrainXY*(i==EL);
    avgStrainYY = avgStrainYY + aStrainYY*(i==EL);
    
    
    
    
end

avgStressXX = avgStressXX';
avgStressXY = avgStressXY';
avgStressYY = avgStressYY';

avgStrainXX = avgStrainXX';
avgStrainXY = avgStrainXY';
avgStrainYY = avgStrainYY';



end