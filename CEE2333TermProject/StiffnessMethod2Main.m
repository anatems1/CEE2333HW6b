%%
% This code is an adaptation of the code from
% http://www.ce.memphis.edu/7117/FEM_codes/Camp's_code/Camp_Matlab_code.html
%
%% Sample FEM code for CEE 2333 Sec 1020
%
%Clear variables and MATLAB workspace
clear all;
clc;
format shortEng;

%% Pre-processing step 
%
% Read inputs from external file
% Input data from Excel spreadsheets 
%
% Input general information

disp("Reading input data");
delimiterIn = ',';
%
%General inputs about the problem
%

% OD = input('Outer radius of the circle (in): ');
% ID = input('Inner radius of the circle (in): ');
OD = 5;
ID = 2;
%ask user for number of elements

% rad_el = input('Number of radial elements: ');
% tang_el = input('Number of tangential elements: ');
rad_el = 5;


%%%%want to make this EVEN numbers onlyy
tang_el = 20;




%ask user for angle of loading

% qmag = input('Magnitude of surface load q(lbs): ');
% qang = input('Input angle (in radians) for surface load q: ');

%ask user for material properities

% E = input('What is the elastic modulus (E - psi) of the material: ');
% v = input('What is the materials poisson ratio (v): ');

NumElem = rad_el*tang_el;
NumDof = NumElem*2;
%determine x and y coordinates of the circle points

xx = zeros(rad_el+1,tang_el);
yy = zeros(rad_el+1,tang_el);

del_rad = 2*pi()/tang_el;
numnodes = (rad_el+1)*(tang_el);
globnodecords = zeros(numnodes,3);
nodenum = 0;

for qz = 1:rad_el+1
    rad1 = OD - ((qz-1)*(OD-ID)/rad_el);
    for jz = 1:tang_el
        ang1 = del_rad*(jz-1);
        xx(qz,jz) = rad1*cos(ang1);
        yy(qz,jz) = rad1*sin(ang1);
        nodenum = nodenum + 1;
        globnodecords(nodenum, 1) = nodenum;
        globnodecords(nodenum, 2) = xx(qz,jz);
        globnodecords(nodenum, 3) = yy(qz,jz);
    end
end

pos1 = 1;
C = zeros(NumElem,4);
for kl = 1:NumElem
        C(kl,1) = pos1;
        C(kl,2) = pos1 + 1;
        C(kl,3) = pos1 + 21;
        C(kl,4) = pos1 + 20;
        pos1 = pos1 + 1;
end
nodenum = nodenum;




%
%% Processing step
%
%Define global stiffness matrix by calling the corresponding stiffness
%function based on element type.
%
fprintf("Forming stiffness matrix\n\n");
fprintf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");
        
    K = StiffnessSpring2(E,C,xx,yy,t,v,NumElem,NumDof);

disp("Applying boundary conditions");
%
%Eliminate rows and columns to reduce the problem to only free DOFs
%
locFree     = find(bcType==0); %Find DOFs that are free, represented by a value of 0 in boundary. See 'Boundary' sheet in the input file.
locFIX  = find(bcType==1); %Find DOFs that are fixed, corresponding to a value of 1. See comment above.
DofFree = length(locFree); %Number of free DOFs

Kelim = zeros(DofFree); %Pre-allocate Kelim, the stiffness matrix corresponding to free DOFs. The function zeros(n) creates an nxn matrix of zeros.
Felim = zeros(DofFree,1); %Pre-allocate Felim, the force vector corresponding to free DOFs. The function zeros(n,m) creates an nxm matrix of zeros.
%
%Copy values of K corresponding to the free DOFs into Kelim.
%Note that the nested loop below could've been replaced with Kelim = K(loc,loc), which would also be faster.
%
fprintf("\n\n The Stiffness matrix after the bc are applied:\n\n");
for i = 1:DofFree %Loop through rows
    for j=1:DofFree %Loop through columns
        Kelim(i,j) = K(locFree(i),locFree(j));  
    end
end

disp(Kelim);
%
%Copy values of F corresponding to the free DOFs into Felim.
%Note that the loop below could've been replaced with Felim(:,1) = F(loc),which would also be faster.
%
for i = 1:DofFree
    Felim(i,1) = bcValue(locFree(i));
end
%
% Solve for free DOFs and store results in d. Note the backslash operator
% will find the fastest way to compute the inverse of Kelim. Typically, it
% will be much faster than calling inv(Kelim) directly.
disp("Finding unknown displacements");
d = Kelim\Felim;
%
%Write full vector of displacements
%
displacement = zeros(NumDof,1); %Pre-allocated a zero-length vector for displacements.
displacement(locFree,1) = d; %Write displacements of free DOFs into the displacement vector - the fixed DOFs have zero displacement.
%

%% Post-processing step
% Evaluate element stresses by calling the appropriate stress function
% corresponding to the element type.

  Se =StressSpring2(E,C,xx,yy,t,v,displacement,NumElem);


%==============================================================

% Write displacements to Excel file
D = zeros(NumNod,3); %Pre-allocate a matrix D to store the displacements. The number of rows is the number of nodes, and there are three columns, as defined below.
%D(i,:) = [ith node number,ith node x-displacement,ith node y-displacement]
for ip = 1:NumNod %Loop through each node
    node = NumDofPerNode*ip; %'node' variable is basically the y-direction DOF at the node ip, this is a not-so-elegant way of expressing DOFs in terms of number of nodes.
    D(ip,1) = ip; %Column 1: Node number
    D(ip,2) = displacement(ip*2-1,1); %Column 2: x-displacement
    D(ip,3) = displacement(ip*2,1); %Column 3: y-displacement
end

disp("Writing the displacements in the Excel file");
%
% 
%
header = {'Displacements'};

% output results
fileName = fullfile(pwd, 'ResultsStiffnessMethod1.csv');
fid = fopen(fileName, 'wt');
% Write headers
fprintf(fid, 'Displacements\n');
fprintf(fid,'%s, %s, %s','Node','u','v');

% Write data.
for i=1:NumNod
    fprintf(fid, '\n%f , %f, %f',D(i,1:3));
end

elemname = zeros(1,NumElem*4);
locnode = zeros(1,NumElem*4);
globnode = zeros(1,NumElem*4);

for jz = 1:NumElem*4 
        elemname(1,jz) = "1";
        locnode(1,jz) = jz;
        globnode(1,jz) = C(1,jz);
end

fprintf(fid, '\nStress\n');
%fprintf(fid,elemname(1:NumElem*4));
fprintf(fid,'%s, %s, %s, %s, %s, %s', 'Element','Gauss Node','Global Node','Stress-x','Stress-y','Stress-xy');

% Write data.
for i=1:NumElem*4
    fprintf(fid,'\n%i ,%i, %i, %f, %f, %f',elemname(1,i),locnode(1,i),globnode(1,i),Se(1,i),Se(2,i),Se(3,i));
 end

fclose(fid);

 