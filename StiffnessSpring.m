function [ K ] = StiffnessSpring( E,C,coorx,coory,area,v,NumEle,NumDof )
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

 %E matrix for plane stress
EMAT = zeros(3,3);
%
    EMAT(1,1) = E;
    EMAT(1,2) = v*E;
    EMAT(1,3) = 0;
	
	EMAT(2,1) = v*E;
    EMAT(2,2) = E;
    EMAT(2,3) = 0;
	
	EMAT(1,1) = 0;
    EMAT(1,2) = 0;
    EMAT(1,3) = (1-v)*E/2;
%HOOKS LAW 2D
EMAT=EMAT*(1/(1-v^2);

%Developing B matrix
B = zeros(3,8);

a = (coorx(1,2) - coorx(1,1))/2;
b = (coory(1,3) - coory(1,2))/2;

for ji = 1:(NumDof/2)
    xp = coorx(1,ji)- coorx(1,1) - a; %values of x' at each node
    yp = coory(1,ji)-coory(1,1) - b; %values of y' at each node
    
    b1 = yp - b;
    b2 = yp + b;
    
    b3 = xp - a;
    b4 = xp + a;
    
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
    B(3,8) = -b4;
    
 
    ans11 = ['x prime NO:',num2str(ji),' = ',num2str(xp)];
    ans12 = ['y prime NO:',num2str(ji),' = ',num2str(yp)];
    ans1 = ['This is the B matrix for Node NO: ',num2str(ji)];
    disp(ans11);
    disp(ans12);
    disp(ans1);
    B = (1/(4*a*b))*B;
    
    disp(B);
end


    




for ii = 1:NumEle %For each element
    
%
%  Determine global coordinates of the nodes of the elements
%

    
    %X3 = coor(C(ii,3)); %x and y coordinates of the third node of the element
	%X4 = coor(C(ii,4)); %x and y coordinates of the fourth node of the element
%
%   Determine the length of the element and the spring stiffness
%
    length = sqrt((X2-X1)^2);
    sk=E(ii)*area(ii)/length;
%
%   Form stiffness matrix in the local coordinates
%
    kk1=[sk,-sk;-sk,sk]; 
  
   
%
    Ktemp = zeros(NumDof);
    for i = 1:2 %Loop through rows corresponding to each DOF
        for j = 1:2 %Loop through columns corresponding to each DOF
           
           %For the ith row
           if i==1 %First node, x-DOF
               index1=C(ii,1);
           elseif i==2 %Second node, x-DOF
               index1=C(ii,2);
           end
           
           %Similarly for the jth column
           if j==1 
               index2=C(ii,1);
           elseif j==2
               index2=C(ii,2);
           end
           
            Ktemp(index1,index2)= kk1(i,j); %Get value in Ke at (i,j) position, and store it in Ktemp at DOFs indicated by (index1,index2). All other values of Ktemp are zero.
        end
    end
    
    %At this point, Ktemp has all the values of Ke (for each element ii)
    %permuted to the appropriate global DOFs. So we can simply add Ktemp
    %into K to assemble the global stiffness matrix.
    K = Ktemp + K;
     

end

