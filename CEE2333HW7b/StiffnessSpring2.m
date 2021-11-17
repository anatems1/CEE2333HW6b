function [ K ] = StiffnessSpring(Emod,C,coorx,coory,area,v,NumEle,NumDof )
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
% K is the global stiffness matrix

K = zeros(NumDof); %Pre-allocated global stiffness matrix of size NumDof x NumDof

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

K = zeros(12,12);
Ktemp = zeros(8,8);
Ktemp_glob =- zeros(12,12);

%transform coordinates
for  i = 1:NumEle
    nodenum = 0;
    locx = zeros(1,4);
    locy = zeros(1,4);
    for ji = 1:4
        locx(1,ji) = coorx(1,C(i,ji));
        locy(1,ji) = coory(1,C(i,ji));
    end
    
    for zt = 1:2
        xi = (t_gs(1,zt))%*2 + 2)/2;
        for xt = 1:2
            nu = (t_gs(1,xt))%*2 + 2)/2;

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

            Ktemp = Wi(1,zt)*Wj(1,zt)*(transpose(B)*EMAT)*B*det(J);
           
            
            %convert local stiffness to global stiffness
            globloc = zeros(1,8);
            for t = 1:4
               globlocx = C(i,t)*2-1;
               globlocy = C(i,t)*2;
               globloc(1,t*2-1) = globlocx;
               globloc(1,t*2) = globlocy;
              
            end   
            %disp(globloc);
             

            Ktemp_glob = zeros(12,12);
            for ti = 1:8 %row
               for tz = 1:8 %column
                    K(globloc(1,ti),globloc(1,tz)) = K(globloc(1,ti),globloc(1,tz)) + Ktemp(ti,tz);
                    Ktemp_glob(globloc(1,ti),globloc(1,tz)) = Ktemp(ti,tz);
               end
            end
            
            nodenum = nodenum + 1;
            %disp(K);
            ans11 = ['nu',num2str(xt),' = ',num2str(nu)];
            ans12 = ['xi',num2str(zt),' = ',num2str(xi)];
            ans3 = ['This is the element B matrix for element NO:',num2str(i),' and gauss Node NO: ',num2str(nodenum)];
            ans1 = ['This is the element stiffness matrix for element NO:',num2str(i),' and gauss Node NO: ',num2str(nodenum)];
            disp(ans11);
            disp(ans12);

            disp(ans3);
            disp(B);
            
            %disp(ans1);
            %disp(Ktemp);
            
            disp(ans1);
            disp(Ktemp_glob);
            
        end
    end

   
end

disp(K);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %Developing B matrix
% B = zeros(3,8);
% 
% a = (coorx(1,2) - coorx(1,1))/2;
% b = (coory(1,3) - coory(1,2))/2;
% 
% 
% %constructing global stiffness matrix
% Wi = [1, 1];
% Wj = [1, 1];
% t_gs = [-1/sqrt(3),1/sqrt(3)];
% 
% curr_sum = 0;
% tot_sum =0;
% nodenum = 1;
% 
% for ji = 1:2
%     %xp = coorx(1,ji)- coorx(1,1) - a; %values of x' at each node
%     xo = (t_gs(1,ji)*2*a + 2*a)/2; %gauss x prime value
%     x_gs = xo - coorx(1,1) - a;
%     
%     Wi1 = Wi(1,ji); %pulling weight value for gauss int.. happens to be 1
%     
%    
%     
%     for jj = 1:2;
%         Wj1 = Wj(1,jj);
%       
%         %yp = coory(1,ji)-coory(1,1) - b; %values of y' at each node
%         yo = (t_gs(1,jj)*2*b + 2*b)/2; %gauss y prime value
%         y_gs = yo - coory(1,1) - b;
%         
%         
%         b1 = y_gs - b;
%         b2 = y_gs + b;
% 
%         b3 = x_gs - a;
%         b4 = x_gs + a;
% 
%         %populate B matrix
%         B(1,1) = b1;
%         B(1,3) = -b1;
%         B(1,5) = b2;
%         B(1,7) = -b2;
% 
%         B(2,2) = b3;
%         B(2,4) = -b4;
%         B(2,6) = b4;
%         B(2,8) = -b3;
% 
%         B(3,1) = b3;
%         B(3,2) = b1;
%         B(3,3) = -b4;
%         B(3,4) = -b1;
%         B(3,5) = b4;
%         B(3,6) = b2;
%         B(3,7) = -b3;
%         B(3,8) = -b2;
%         
%         ans11 = ['x prime NO:',num2str(ji),' = ',num2str(x_gs)];
%         ans12 = ['y prime NO:',num2str(jj),' = ',num2str(y_gs)];
%         ans1 = ['This is the B matrix for gauss Node NO: ',num2str(nodenum)];
%         disp(ans11);
%         disp(ans12);
%         disp(ans1);
%         
%         
%         B = (1/(4*a*b))*B;
%         disp(B);
%         
%         
%         curr_sum = 1*a*b*Wi1*Wj1*transpose(B)*EMAT*B;%
%         
%         ans2 = ['The global matrix for gauss Node NO: ',num2str(nodenum)];
%         disp(ans2);
%         disp(curr_sum);
%         
%         nodenum = nodenum + 1;
%         
%         tot_sum = tot_sum + curr_sum;%
%     end
%  
%     
% 
% end
% 
% K = tot_sum;







end

