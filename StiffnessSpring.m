function [ K ] = StiffnessSpring( E,C,coor,area,v,NumEle,NumDof )
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

for ii = 1:NumEle %For each element
    
%
%  Determine global coordinates of the nodes of the elements
%
    X1 = coor(C(ii,1)); %x and y coordinates of the first node of the element
    X2 = coor(C(ii,2)); %x and y coordinates of the second node of the element
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

