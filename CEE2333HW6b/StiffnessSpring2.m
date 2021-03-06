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

%Developing B matrix
B = zeros(3,8);

a = (coorx(1,2) - coorx(1,1))/2;
b = (coory(1,3) - coory(1,2))/2;


%constructing global stiffness matrix
Wi = [1, 1];
Wj = [1, 1];
t_gs = [-1/sqrt(3),1/sqrt(3)];

curr_sum = 0;
tot_sum =0;
nodenum = 1;

for ji = 1:2
    %xp = coorx(1,ji)- coorx(1,1) - a; %values of x' at each node
    xo = (t_gs(1,ji)*2*a + 2*a)/2; %gauss x prime value
    x_gs = xo - coorx(1,1) - a;
    
    Wi1 = Wi(1,ji); %pulling weight value for gauss int.. happens to be 1
    
   
    
    for jj = 1:2;
        Wj1 = Wj(1,jj);
      
        %yp = coory(1,ji)-coory(1,1) - b; %values of y' at each node
        yo = (t_gs(1,jj)*2*b + 2*b)/2; %gauss y prime value
        y_gs = yo - coory(1,1) - b;
        
        
        b1 = y_gs - b;
        b2 = y_gs + b;

        b3 = x_gs - a;
        b4 = x_gs + a;

        %populate B matrix
        B(1,1) = b1;
        B(1,3) = -b1;
        B(1,5) = b2;
        B(1,7) = -b2;

        B(2,2) = b3;
        B(2,4) = -b4;
        B(2,6) = b4;
        B(2,8) = -b3;

        B(3,1) = b3;
        B(3,2) = b1;
        B(3,3) = -b4;
        B(3,4) = -b1;
        B(3,5) = b4;
        B(3,6) = b2;
        B(3,7) = -b3;
        B(3,8) = -b2;
        
        ans11 = ['x prime NO:',num2str(ji),' = ',num2str(x_gs)];
        ans12 = ['y prime NO:',num2str(jj),' = ',num2str(y_gs)];
        ans1 = ['This is the B matrix for gauss Node NO: ',num2str(nodenum)];
        disp(ans11);
        disp(ans12);
        disp(ans1);
        
        
        B = (1/(4*a*b))*B;
        disp(B);
        
        
        curr_sum = 1*a*b*Wi1*Wj1*transpose(B)*EMAT*B;%
        
        ans2 = ['The global matrix for gauss Node NO: ',num2str(nodenum)];
        disp(ans2);
        disp(curr_sum);
        
        nodenum = nodenum + 1;
        
        tot_sum = tot_sum + curr_sum;%
    end
 
    

end

K = tot_sum;





% for ii = 1:NumEle %For each element
%     
% %
% %  Determine global coordinates of the nodes of the elements
% %
% 
%     
%     %X3 = coor(C(ii,3)); %x and y coordinates of the third node of the element
% 	%X4 = coor(C(ii,4)); %x and y coordinates of the fourth node of the element
% %
% %   Determine the length of the element and the spring stiffness
% %
%     length = sqrt((X2-X1)^2);
%     sk=Emod(ii)*area(ii)/length;
% %
% %   Form stiffness matrix in the local coordinates
% %
%     kk1=[sk,-sk;-sk,sk]; 
%   
%    
% %
%     Ktemp = zeros(NumDof);
%     for i = 1:2 %Loop through rows corresponding to each DOF
%         for j = 1:2 %Loop through columns corresponding to each DOF
%            
%            %For the ith row
%            if i==1 %First node, x-DOF
%                index1=C(ii,1);
%            elseif i==2 %Second node, x-DOF
%                index1=C(ii,2);
%            end
%            
%            %Similarly for the jth column
%            if j==1 
%                index2=C(ii,1);
%            elseif j==2
%                index2=C(ii,2);
%            end
%            
%             Ktemp(index1,index2)= kk1(i,j); %Get value in Ke at (i,j) position, and store it in Ktemp at DOFs indicated by (index1,index2). All other values of Ktemp are zero.
%         end
%     end
%     
%     %At this point, Ktemp has all the values of Ke (for each element ii)
%     %permuted to the appropriate global DOFs. So we can simply add Ktemp
%     %into K to assemble the global stiffness matrix.
%     K = Ktemp + K;
     

end

