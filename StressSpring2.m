% function [ Se ] = StressSpring( E,C,coor,area,v,displacement,NumEle )
% % Stress for each rod element [Se]=[D][B]{d}

% % Inputs
% % E is the modulus of elastcity of the element
% % C is the connectivity matrix
% % s is the matrix of coordinates of the nodes
% % v is Poisson's ratio
% % displacement is the global displacement vector (calculated from the FEM
% % problem)
% % NumEle is the number of elements

% % Outputs
% % Se is the element stress matrix, each row consists of three entries: x
% % stress as well as zero values for y stress and xy stress for that
% % element.

% Se = zeros(NumEle,3); %Pre-allocated Se array, as defined above
% A = zeros(NumEle,1);
% fg = zeros(2,1);
% uu1= zeros(2,1);

% for ii = 1:NumEle %For each element
% %
% %  Determine global coordinates of the nodes of the elements
% %
    % S1 = coor(C(ii,1)); %x and y coordinates of the first node of the element
    % S2 = coor(C(ii,2)); %x and y coordinates of the second node of the element
% %
% %   Determine the length of the element and the spring stiffness
% %
    % length = sqrt((S2-S1)^2);
    % sk=E(ii)*area(ii)/length;
% %
% %   Form stiffness matrix in the local coordinates
% %
    % kk1=[sk,-sk;-sk,sk]; 
   % dtemp(1:2) = 0;              
   % for i = 1:2  %Loop through each DOF of the element
       
           % %For the ith DOF
           % if i==1 %First node, x-DOF (global DOF)
               % index1=C(ii,1);
           % elseif i==2 %First node, y-DOF (global DOF)
               % index1=C(ii,2);
           % end
           
            % dtemp(i)= displacement(index1,1); %Get displacement at global DOF index1 and store it in local DOF i
            
   % end
 % %
 % %  Compute the element nodal forces in the global coordinates
 % %
    % fg = kk1*transpose(dtemp);
 % %
 % %  Compute element stress
 % %
    % Se(ii,1) =fg(2)/area(ii) ; %Evaluate stress for element ii
     
% end

