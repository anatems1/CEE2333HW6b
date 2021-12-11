function [ bcValue ] = circloading(Emod,C,coorx,coory,area,v,NumEle,NumDof,qmag,qang,bcValue,partial_load,xa,ya,choice,tang_el,lock)
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

%test
% Outputs
%% 
% K is the global stiffness matrix

bcValue = zeros(NumDof,1); %Pre-allocated global stiffness matrix of size NumDof x NumDof
bcTemp = zeros(4,1);

Wi = [1, 1];
Wj = [1, 1];
t_gs = [-1/sqrt(3),1/sqrt(3)];



%transform coordinates
for  i = 1:NumEle
    bda = 0;
    bdb = 0;
    nodenum = 0;
    locx = zeros(1,4);
    locy = zeros(1,4);
    
    for ji = 1:4
        locx(1,ji) = coorx(1,C(i,ji));
        locy(1,ji) = coory(1,C(i,ji));
    end
    
xa_loc = 0;

%determine upper and lower bounds for gauss int along surface
for ii = 1:4
    if(i == partial_load(1,ii))
        x1 =  locx(1,1) + locx(1,2) + locx(1,3) + locx(1,4);
        xb = (-locx(1,1) - locx(1,2) + locx(1,3) + locx(1,4));
        xc = (-locx(1,1) + locx(1,2) + locx(1,3) - locx(1,4));
        xd = (locx(1,1) - locx(1,2) + locx(1,3) - locx(1,4));
        x1 = x1 * 0.25;
        xb = xb * 0.25;
        xc = xc * 0.25;
        xd = xd * 0.25;

        xi_pt = (xa(1,ii)-x1-(xb*-1))/(xc+(xd*-1));
        xa_loc = xa(1,ii);
        if lock == 1 
            bda = -1;
            bdb = 1;
        else
            if xi_pt < 0 
                bdb = xi_pt;
                bda = -1;

            elseif xi_pt>=0 
                bdb = 1;
                bda = xi_pt; 
            end
        end
    elseif (i > partial_load(1,1) && i < partial_load(1,2)) || (i > partial_load(1,3) && i < partial_load(1,4))
        bda = -1;
        bdb = 1;
    elseif (choice == 2 && i > partial_load(1,1))
        bda = -1;
        bdb = 1;
    end
end


%determine length of member (surface)
%L = sqrt((coorx(1,C(i,2)) - coorx(1,C(i,1)))^2 + (coory(1,C(i,2)) - coory(1,C(i,1)))^2);
L = abs(coorx(1,C(i,2)) - coorx(1,C(i,1)));

q = 0;

%determine direction of loading
if(i >= partial_load(1,1) && i <= partial_load(1,2) && choice == 1)
    q = -qmag; %adjust to angle
        
elseif(i >= partial_load(1,3) && i <= partial_load(1,4) && choice == 1)
    q = qmag; %adjust to angle

elseif choice == 2 && i >= partial_load(1,1) && i < tang_el/4 +1
    q = -qmag;
    
else
    q = 0;

end

%intergrate and determine loading at nodes
    %xi = ((t_gs(1,xt))*(bdb-bda) + (bdb-bda))/2;

    %N1 = 0.25*(1-xi)*(1-(-1));
    %N2 = 0.25*(1+xi)*(1-(-1));
    int_N1 = 0.5*((bdb - bdb^2/2) - (bda - bda^2/2)); %integrated shape functions
    int_N2 = 0.5*((bdb + bdb^2/2) - (bda + bda^2/2));

    N = [int_N1, 0, int_N2, 0 ; 0, int_N1, 0, int_N2]; %surface shape function matrix
    N = transpose(N);

    q_mat = [0; q];

    bcTemp = 1*L/2*N*q_mat;    

    bcValue(C(i,1)*2,1) = bcValue(C(i,1)*2,1) + bcTemp(2,1);
    if choice == 2 && i == tang_el/4
        bcValue(C(i,2)*2,1) = bcValue(C(i,2)*2,1) + bcTemp(4,1)*2;       
    else
        bcValue(C(i,2)*2,1) = bcValue(C(i,2)*2,1) + bcTemp(4,1);
    end 
end

if choice ==2 
bcValue(2*(tang_el/4 + 1),1) = 0.5*bcValue(2*(tang_el/4 + 1),1);
end
bcValue = bcValue;




