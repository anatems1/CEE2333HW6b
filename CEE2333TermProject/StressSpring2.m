function [ Se ] = StressSpring( Emod,C,coorx,coory,area,v,displacement,NumEle )
% Stress for each rod element [Se]=[D][B]{d}

% Inputs
% E is the modulus of elastcity of the element
% C is the connectivity matrix
% s is the matrix of coordinates of the nodes
% v is Poisson's ratio
% displacement is the global displacement vector (calculated from the FEM
% problem)
% NumEle is the number of elements

% Outputs
% Se is the element stress matrix, each row consists of three entries: x
% stress as well as zero values for y stress and xy stress for that
% element.

% Se = zeros(NumEle,3); %Pre-allocated Se array, as defined above

Se1 = zeros(3,NumEle*2);

% A = zeros(NumEle,1);
% fg = zeros(2,1);
% uu1= zeros(2,1);

Emod = Emod(1,1);
vpois = v(1,1);
%E matrix for plane stress
EMAT = zeros(3,3);
%
    EMAT(1,1) = Emod;
    EMAT(1,2) = vpois*Emod;
    EMAT(1,3) = 0;
	
	EMAT(2,1) = vpois*Emod;
    EMAT(2,2) = Emod;
    EMAT(2,3) = 0;
	
	EMAT(3,1) = 0;
    EMAT(3,2) = 0;
    EMAT(3,3) = (1-vpois)*Emod/2;
%HOOKS LAW 2D
EMAT = EMAT*(1/(1-vpois^2));

fprintf("\nThe Plane Stress [E] matrix is:\n");
disp(EMAT);

Wi = [1, 1];
Wj = [1, 1];
t_gs = [-1/sqrt(3),1/sqrt(3)];

nodenum = 0;

%transform coordinates
for  i = 1:NumEle
    locx = zeros(1,4);
    locy = zeros(1,4);
    for ji = 1:4
        locx(1,ji) = coorx(1,C(i,ji));
        locy(1,ji) = coory(1,C(i,ji));
    end
    
    for zt = 1:2
        nu = (t_gs(1,zt));%*2 + 2)/2;
        for xt = 1:2
            xi = (t_gs(1,xt));%*2 + 2)/2;

            x = locx(1,1)*0.25*(1-xi)*(1-nu) + locx(1,2)*0.25*(1+xi)*(1-nu) + locx(1,3)*0.25*(1+xi)*(1+nu) + locx(1,4)*0.25*(1-xi)*(1+nu);
            y = locy(1,1)*0.25*(1-xi)*(1-nu) + locy(1,2)*0.25*(1+xi)*(1-nu) + locy(1,3)*0.25*(1+xi)*(1+nu) + locy(1,4)*0.25*(1-xi)*(1+nu);

            %del_x/del_xi = (a + bnu + cxi + dnu*xi)
            xb = (-locx(1,1) - locx(1,2) + locx(1,3) + locx(1,4));
            xc = (-locx(1,1) + locx(1,2) + locx(1,3) - locx(1,4));
            xd = (locx(1,1) - locx(1,2) + locx(1,3) - locx(1,4));
            xb = xb * 0.25;
            xc = xc * 0.25;
            xd = xd * 0.25;

            yb = (-locy(1,1) - locy(1,2) + locy(1,3) + locy(1,4));
            yc = (-locy(1,1) + locy(1,2) + locy(1,3) - locy(1,4));
            yd = (locy(1,1) - locy(1,2) + locy(1,3) - locy(1,4));
            yb = yb * 0.25;
            yc = yc * 0.25;
            yd = yd * 0.25;

            J = zeros(2,2);
            J(1,1) = xc + xd*nu;
            J(1,2) = yc + yd*nu;
            J(2,1) = xb + xd*xi;
            J(2,2) = yb + yd*xi;

            gam = zeros(2,2);
            gam = inv(J);

            %developing B matrix
            B1 = zeros(3,4);
            B1(1,1) = 1;
            B1(3,2) = 1;
            B1(3,3) = 1;
            B1(2,4) = 1;

            B2 = zeros(4,4);
            B2(1:2,1:2) = gam;
            B2(3:4,3:4) = gam;

            b31 = -1 + nu;
            b32 = 1 + nu;
            b33 = -1 + xi;
            b34 = 1 + xi;

            B3 = zeros(3,8);
            B3(1,1) = b31;
            B3(1,3) = -b31;
            B3(1,5) = b32;
            B3(1,7) = -b32;
            B3(2,1) = b33;
            B3(2,3) = -b34;
            B3(2,5) = b34;
            B3(2,7) = -b33;

            B3(3:4,2:8) = B3(1:2,1:7);
            B3 = 0.25*B3;

            B = B1*B2*B3;

            %convert local stiffness to global stiffness
            globloc = zeros(1,8);
            for t = 1:4
               globlocx = C(i,t)*2-1;
               globlocy = C(i,t)*2;
               globloc(1,t*2-1) = globlocx;
               globloc(1,t*2) = globlocy;
              
            end   
            

%the only difference in code from the stiffnesssping code below this
%line\/\//\/\//\/\/\/
            nodenum = nodenum + 1;
            disp1 = zeros(8,1);
            %pull displacements out for element
            for tj = 1:8
                disp1(tj,1) = displacement(globloc(1,tj),1);
            end
            
            Se1(1:3,nodenum) = EMAT * B * disp1;
            
        end
    end
end

Se = Se1;

%disp(Se);

%disp('test');
% %
% %  Determine global coordinates of the nodes of the elements
% %
%     S1 = coor(C(ii,1)); %x and y coordinates of the first node of the element
%     S2 = coor(C(ii,2)); %x and y coordinates of the second node of the element
% %
% %   Determine the length of the element and the spring stiffness
% %
%     length = sqrt((S2-S1)^2);
%     sk=E(ii)*area(ii)/length;
% %
% %   Form stiffness matrix in the local coordinates
% %
%     kk1=[sk,-sk;-sk,sk]; 
%    dtemp(1:2) = 0;              
%    for i = 1:2  %Loop through each DOF of the element
%        
%            %For the ith DOF
%            if i==1 %First node, x-DOF (global DOF)
%                index1=C(ii,1);
%            elseif i==2 %First node, y-DOF (global DOF)
%                index1=C(ii,2);
%            end
%            
%             dtemp(i)= displacement(index1,1); %Get displacement at global DOF index1 and store it in local DOF i
%             
%    end
%  %
%  %  Compute the element nodal forces in the global coordinates
%  %
%     fg = kk1*transpose(dtemp);
%  %
%  %  Compute element stress
%  %
%     Se(ii,1) =fg(2)/area(ii) ; %Evaluate stress for element ii