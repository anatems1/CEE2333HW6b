function [ bcValue ] = circloading(Emod,C,coorx,coory,area,v,NumEle,NumDof,qmag,qang,bcValue,partial_load,xa,ya)
%
% Stiffness matrix for each rod element 
%
% Inputs
% E is the modulus of elastcity of the element
% C is the connectivity matrix
% coor is the matrix of coordinates of the nodes
% area is the cross-section area of each element
% v is Poisson's ratio
% NumEle is the number of elements
% NumDof is the number of DOFs

% Outputs
%% 
% K is the global stiffness matrix

bcValue = zeros(NumDof); %Pre-allocated global stiffness matrix of size NumDof x NumDof


Wi = [1, 1];
Wj = [1, 1];
t_gs = [-1/sqrt(3),1/sqrt(3)];

bda = 0;
bdb = 0;

%transform coordinates
for  i = 1:NumEle
    nodenum = 0;
    locx = zeros(1,4);
    locy = zeros(1,4);
    
    for ji = 1:4
        locx(1,ji) = coorx(1,C(i,ji));
        locy(1,ji) = coory(1,C(i,ji));
    end
    
%del_x/del_xi = (a + bnu + cxi + dnu*xi)
x1 =  locx(1,1) + locx(1,2) + locx(1,3) + locx(1,4);
xb = (-locx(1,1) - locx(1,2) + locx(1,3) + locx(1,4));
xc = (-locx(1,1) + locx(1,2) + locx(1,3) - locx(1,4));
xd = (locx(1,1) - locx(1,2) + locx(1,3) - locx(1,4));
x1 = x1 * 0.25;
xb = xb * 0.25;
xc = xc * 0.25;
xd = xd * 0.25;

xi1 = 
    
%determine upper and lower bounds for gauss int along surface
for ii = 1:4
    if(i = partial_load(1,i))
        bda = 
        bdb = 
    else
        bda = 1;
        bdb = -1;
    end
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
            
            %convert local stiffness to global stiffness
            globloc = zeros(1,8);
            for t = 1:4
               globlocx = C(i,t)*2-1;
               globlocy = C(i,t)*2;
               globloc(1,t*2-1) = globlocx;
               globloc(1,t*2) = globlocy;
              
            end   
            %disp(globloc);
            
        end
    end 
end

fprintf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");

fprintf("The Global Stiffness Matrix is: \n\n");
%% determine which output we are selecting
    K = K

end





