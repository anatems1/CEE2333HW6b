%%
% This code is an adaptation of the code from
% http://www.ce.memphis.edu/7117/FEM_codes/Camp's_code/Camp_Matlab_code.html
%
%% Sample FEM code for CEE 2333 Sec 1020
%
%Clear variables and MATLAB workspace
clear all;
clc;
close all;
format shortEng;

%test
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

%inputs used during the creation of this program
% OD = 6;
% ID = 3;
% qmag = 1000;
% qang = 45;
% E = 36000;
% v=0.15;
% rad_el = 5;
% tang_el = 20;
% defscale = 1;


OD = input('Outer radius of the circle (in): ');
ID = input('Inner radius of the circle (in): ');

% %ask user for angle and magnitude of loading
qmag = input('Magnitude of surface load q(lbs): ');
qang = input('Input angle (in radians) for surface load q: ');

ang_check = 0; %checks that the angle of loading is not exceeding the bounds of half of the total surface
while ang_check < 0
    if qang > 180 || qang < 0
        fprintf("\nPlease enter an angle between 0 and 180 \n");
        qang = input('Input angle (in radians) for surface load q: ');
    else
        ang_check = 1;
    end
end

%ask user for material properities
E = input('What is the elastic modulus (E - psi) of the material: ');
v = input('What is the materials poisson ratio (v): ');

%ask user for number of elements
rad_el = input('Number of radial elements: ');
tang_el = input('Number of tangential elements: ');

%%%%want to make this EVEN numbers onlyy and the number must be greater than 4
ev_check = 0;
while ev_check < 1
    if mod(tang_el,2) > 0
        fprintf("\nPlease enter an even number of elements!!! \n");
        tang_el = input('Number of tangential elements: ');
    elseif tang_el <= 4 && choice == 1
        fprintf("\nPlease enter a number greater than or equal to 4!!! \n");
        tang_el = input('Number of tangential elements: ');
    elseif mod(tang_el,4)>0 && choice == 2
        fprintf("\nSince the analysis is for a quarter circle, the number of elements must be in multiples of 4.... \n");
        tang_el = input('Number of tangential elements: ');
    else
        ev_check = 2;
    end
end

%asks user for a deformation scale that will be used in the postprocessor
defscale = input('Enter deformation scale: ');

qang = qang * pi()/180;%converts angle to radians
t=1;%the elements is always a thickness of 1 for plane stress

%%asks user if the analysis should be a quarter or the full circle
choice = menu('Please choose analysis option:', 'Full Circle', 'Quarter Circle');
intchoice = menu('Please choose integration type:', 'Full integration', 'Reduced integration');
postchoice = menu('Would you like to plot stress w/ colormap?(takes much longer to process)','Yes','No');
%determines coordinates for full or quarter circle
switch(choice)
    case 1 %full
        NumElem = rad_el*tang_el;
        numnodes = (rad_el+1)*(tang_el);
        NumDof = numnodes*2;
        loop2 = tang_el;
    case 2 %quarter circle
        NumElem = rad_el*tang_el/4;
        numnodes = (rad_el+1)*((tang_el)/4+1);
        NumDof = numnodes*2;
        loop2 = tang_el/4 + 1;
end


%determine x and y coordinates of the circle points
xx = zeros(rad_el+1,loop2);
yy = zeros(rad_el+1,loop2);


del_rad = 2*pi()/tang_el;

globnodecords = zeros(numnodes,3);
nodenum = 0;

bcType = zeros(NumDof,1); %boundary type --> 1 = fixed, 0 = free (forces)
bcValue = zeros(NumDof,1); %forces



%determines where the bounds of the partial surface loading is occuring
xa = zeros(1,4);
ya = zeros(1,4);
partial_load = zeros(1,4);

xa(1,1) = OD*cos(pi()/2 - qang/2);
ya(1,1) = OD*sin(pi()/2 - qang/2);
xa(1,2) = -1*xa(1,1);
ya(1,2) = ya(1,1);
xa(1,3) = xa(1,2);
ya(1,3) = -1*ya(1,1);
xa(1,4) = xa(1,1);
ya(1,4) = ya(1,3);

temp_el = 0;
loc1 = 1;
lock = 0; %this lock determines wether or no the partial loading will occur, if not the elements will recieve half of the contribution

for qz = 1:rad_el+1
    rad1 = OD - ((qz-1)*(OD-ID)/rad_el);
    for jz = 1:loop2
        ang1 = del_rad*(jz-1);
        xx(qz,jz) = rad1*cos(ang1);
        yy(qz,jz) = rad1*sin(ang1);
        nodenum = nodenum + 1;
        temp_el = temp_el + 1;
        globnodecords(nodenum, 1) = nodenum;;
        globnodecords(nodenum, 2) = xx(qz,jz);
        globnodecords(nodenum, 3) = yy(qz,jz);
        
        %check if angle of end of pressure is in between nodes
        if ang1 == pi()/2 - qang/2 || ang1 == 3*pi()/2 - qang/2  && rad1 == OD;
            partial_load(1,loc1) = temp_el;
            loc1 = loc1 + 1;
            lock = 1; %lock = 1 tells us that there is no partial loading
        elseif ang1 == pi()/2 + qang/2 || ang1 == 3*pi()/2 + qang/2 && rad1 == OD;
            partial_load(1,loc1) = temp_el-1;
            loc1 = loc1 + 1;
            lock = 1;
        elseif lock == 0; %if there is partial loading
            if (ang1 ~= (3*pi()/2 - qang/2)) && ((ang1 > (pi()/2 - qang/2)) && (ang1-del_rad < (pi()/2 - qang/2))) && rad1 == OD;
                partial_load(1,loc1) = temp_el-1;
                loc1 = loc1 + 1;
            elseif (ang1 ~= (3*pi()/2 + qang/2)) && ((ang1 > (pi()/2 + qang/2)) && (ang1-del_rad < (pi()/2 + qang/2))) && rad1 == OD;
                partial_load(1,loc1) = temp_el-1;
                loc1 = loc1 + 1;
            elseif (ang1 ~= (3*pi()/2 - qang/2)) && ((ang1 > (3*pi()/2 - qang/2)) && (ang1-del_rad < (3*pi()/2 - qang/2))) && rad1 == OD;
                partial_load(1,loc1) = temp_el-1;
                loc1 = loc1 + 1;
            elseif (ang1 ~= (3*pi()/2 + qang/2)) && ((ang1 > (3*pi()/2 + qang/2)) && (ang1-del_rad < (3*pi()/2 + qang/2))) && rad1 == OD;
                partial_load(1,loc1) = temp_el-1;
                loc1 = loc1 + 1;
            end
        end
        
        %establishing boundary conditions for the structure
        if (choice == 1) && (xx(qz,jz) == OD || xx(qz,jz) == -OD)
            if xx(qz,jz) == OD
                bcType(nodenum*2) = 1;
            elseif xx(qz,jz) == -OD
                bcType(nodenum*2) = 1;
                bcType(nodenum*2-1) = 1;
            end
        elseif choice == 2
            if ang1 == 0
                bcType(nodenum*2) = 1;
            elseif ang1 == pi()/2
                bcType(nodenum*2-1) = 1;
            end
        end
    end
end

%creating the connectivity array used in further proccesses
pos1 = 1;
C = zeros(NumElem,4);
for kl = 1:NumElem
    C(kl,1) = pos1;
    C(kl,4) = pos1 + loop2;
    if mod(kl,tang_el) == 0 && choice == 1
        C(kl,2) = pos1 + 1 - loop2;
        C(kl,3) = pos1 + 1;
    else
        C(kl,2) = pos1 + 1;
        C(kl,3) = pos1 + loop2 + 1;
    end
    if choice == 1
        pos1 = pos1 + 1;
    elseif choice == 2 && mod(pos1+1,loop2) == 0
        pos1 = pos1 + 2;
    else
        pos1 = pos1 + 1;
    end
end

coorx = zeros(numnodes,1);
coory = zeros(numnodes,1);

coorx(1:numnodes,1) = globnodecords(1:numnodes,2);
coory(1:numnodes,1) = globnodecords(1:numnodes,3);

coorx = transpose(coorx);
coory = transpose(coory);

%
%% Processing step
%
%Define global stiffness matrix by calling the corresponding stiffness
%function based on element type.
%
fprintf("Forming stiffness matrix\n\n");
fprintf("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n");

K = StiffnessSpring2(E,C,coorx,coory,t,v,NumElem,NumDof,intchoice);
bcValue = circloading(E,C,coorx,coory,t,v,tang_el,NumDof,qmag,qang,bcValue,partial_load,xa,ya,choice,tang_el,lock); %to create the force vector with surface loads
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

%disp(Kelim);
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

Se =StressSpring2(E,C,coorx,coory,t,v,displacement,NumElem,intchoice);


%==============================================================

% Write displacements to Excel file
D = zeros(numnodes,3); %Pre-allocate a matrix D to store the displacements. The number of rows is the number of nodes, and there are three columns, as defined below.
%D(i,:) = [ith node number,ith node x-displacement,ith node y-displacement]
for ip = 1:numnodes %Loop through each node
    node = 2*ip; %'node' variable is basically the y-direction DOF at the node ip, this is a not-so-elegant way of expressing DOFs in terms of number of nodes.
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
if choice == 1
    fileName = fullfile(pwd, 'ResultsStiffnessMethod_FullCircle.csv');
else
    fileName = fullfile(pwd, 'ResultsStiffnessMethod_QuarterCircle.csv');
end

fid = fopen(fileName, 'wt');
% Write headers
fprintf(fid, 'Displacements\n');
fprintf(fid,'%s, %s, %s','Node','u','v');

defcoorx = zeros(1,numnodes);
defcoory = zeros(1,numnodes);

% Write data.
for i=1:numnodes
    fprintf(fid, '\n%f , %f, %f',D(i,1:3));
    defcoorx(1,i) = coorx(1,i) + defscale*D(i,2);
    defcoory(1,i) = coory(1,i) + defscale*D(i,3);
end



elemname = zeros(1,NumElem*4);
locnode = zeros(1,NumElem*4);
globnode = zeros(1,NumElem*4);

for jz = 1:NumElem
    for jj = 1:4
        elemname(1,jz*4-4+jj) = jz;
        locnode(1,jz*4-4+jj) = jj;
        globnode(1,jz*4-4+jj) = C(jz,jj);
    end
end

fprintf(fid, '\nStress\n');
%fprintf(fid,elemname(1:NumElem*4));
fprintf(fid,'%s, %s, %s, %s, %s, %s', 'Element','Gauss Node','Global Node','Stress-x','Stress-y','Stress-xy');

% Write data.
for i=1:NumElem*4
    fprintf(fid,'\n%i ,%i, %i, %f, %f, %f',elemname(1,i),locnode(1,i),globnode(1,i),Se(1,i),Se(2,i),Se(3,i));
end

fclose(fid);




%arrays and variables used in the plotting processes below
xz = zeros(1,5);
yz = zeros(1,5);
stress1 = zeros(1,5);

xx1 = zeros(1,2);
yy1 = zeros(1,2);
mx1 = 0;

%%
%plotting the data
if choice == 1
    figure('Name','Full Circle','NumberTitle','off')
    lim = OD*1.25;
    lim2 = -lim;
    mx1 = tang_el;
    dim1 = 50;
else
    figure('Name','Quarter Circle','NumberTitle','off')
    lim = OD*1.25;
    lim2 = -1;
    mx1 = tang_el/4+1;
    dim1 = 50
end

if postchoice == 1 %determines wether or not stress plotting will show
    set(gcf,'position',[0,dim1,1800,800],'color','w')
else
    set(gcf,'position',[0,dim1,1800,400],'color','w')
end


if postchoice == 1%determines wether or not stress plotting will show
    subplot(2,3,1)
else
    subplot(1,3,1)
end

%original model
for z = 1:NumElem
    for ii = 1:4
        xz(1,ii) = coorx(1,C(z,ii));
        yz(1,ii) = coory(1,C(z,ii));
    end
    xz(1,5) = coorx(1,C(z,1));
    yz(1,5) = coory(1,C(z,1));
    
    %plotting each element
    plot(xz,yz,"-o",'color',[0.3010 0.7450 0.9330],'MarkerSize',1,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    hold on
end

labels = zeros(1,numnodes);
for z=1:numnodes
    labels(1,z) = z;
end

%creating node labels
labels = string(labels);

title('Original Model')
plot(coorx,coory,".",'MarkerFaceColor','black','MarkerEdgeColor','black')
text(coorx,coory,labels,'HorizontalAlignment','center',...
    'VerticalAlignment','top','fontsize',6)


xlim([lim2 lim])
ylim([lim2 lim])


%creating loading lines on original plot
for lo = 1:2
    if lo == 1
        ang3 = pi()/2-qang/2;
    else
        ang3 = pi()/2+qang/2;
    end
    
    xx1(1,1) = (OD+1)*cos(ang3);
    xx1(1,2) = -(OD+1)*cos(ang3);
    yy1(1,1) = (OD+1)*sin(ang3);
    yy1(1,2) = -(OD+1)*sin(ang3);
    plot(xx1,yy1,"-",'color',[0.5 0.2 0.5])
end
hold off


%force vectors
if postchoice == 1
    subplot(2,3,2)
else
    subplot(1,3,2)
end

for z = 1:NumElem
    for ii = 1:4
        xz(1,ii) = coorx(1,C(z,ii));
        yz(1,ii) = coory(1,C(z,ii));
    end
    xz(1,5) = coorx(1,C(z,1));
    yz(1,5) = coory(1,C(z,1));
    
    plot(xz,yz,"-o",'color',[0.3010 0.7450 0.9330],'MarkerSize',1,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    hold on
end

title('Force Vectors')
xlim([lim2 lim])
ylim([lim2 lim])



%get location of equivalent nodal forces and plot on graph
load_pts = zeros(1,tang_el);
count = 0;

for  i = 1:mx1
    for ii = 1:4
        if(i == partial_load(1,ii) || i == partial_load(1,ii) + 1 && choice ==1)
            count = count + 1;
            load_pts(1,count) = i;
        end
    end
    if (i > partial_load(1,1)+1 && i < partial_load(1,2)) || (i > partial_load(1,3)+1 && i < partial_load(1,4) && choice ==1)
        count = count + 1;
        load_pts(1,count) = i;
    elseif (choice == 2 && i > partial_load(1,1))
        count = count + 1;
        load_pts(1,count) = i;
    end
    
end

loadx = zeros(1,count);
loady = zeros(1,count);
loadmagy = zeros(1,count);

for i = 1:count
    loadx(1,i) = coorx(1,load_pts(1,i));
    loady(1,i) = coory(1,load_pts(1,i));
    loadmagy(1,i) = round(abs(bcValue(load_pts(1,i)*2,1)));
end
loadmagy = string(loadmagy);

if choice == 1
    plot(loadx(1,1:(count/2)),loady(1,1:count/2),'.')
    text(loadx(1,1:(count/2)),loady(1,1:count/2),'$\downarrow $','Interpreter','latex','HorizontalAlignment','center',...
        'VerticalAlignment','top','fontsize',12,'fontweight','bold')
    text(loadx(1,1:(count/2)),loady(1,1:count/2),loadmagy(1,1:count/2),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',8)
    plot(loadx(1,(count/2+1):count),loady(1,count/2+1:count),'.')
    text(loadx(1,(count/2+1):count),loady(1,count/2+1:count),'$\uparrow $','Interpreter','latex','HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',12,'fontweight','bold')
    text(loadx(1,(count/2+1):count),loady(1,count/2+1:count),loadmagy(1,1:count/2),'HorizontalAlignment','center',...
        'VerticalAlignment','top','fontsize',8)
else
    plot(loadx(1,1:count),loady(1,1:count),'.')
    text(loadx(1,1:count),loady(1,1:count),'$\downarrow $','Interpreter','latex','HorizontalAlignment','center',...
        'VerticalAlignment','top','fontsize',12)
    text(loadx(1,1:count),loady(1,1:count),loadmagy(1,1:count),'HorizontalAlignment','center',...
        'VerticalAlignment','bottom','fontsize',8)
    
end




%displacement field
if postchoice == 1
    subplot(2,3,3)
else
    subplot(1,3,3)
end

%plotting shadow of original strucutre
for z = 1:NumElem
    for ii = 1:4
        xz(1,ii) = coorx(1,C(z,ii));
        yz(1,ii) = coory(1,C(z,ii));
    end
    xz(1,5) = coorx(1,C(z,1));
    yz(1,5) = coory(1,C(z,1));
    
    plot(xz,yz,"-",'linewidth',2,'color',[0.84 0.84 0.84])
    hold on
end

%plotting deformed shape
for z = 1:NumElem
    for ii = 1:4
        xz(1,ii) = defcoorx(1,C(z,ii));
        yz(1,ii) = defcoory(1,C(z,ii));
    end
    xz(1,5) = defcoorx(1,C(z,1));
    yz(1,5) = defcoory(1,C(z,1));
    plot(xz,yz,"-o",'color',[1 0 0],'MarkerSize',3,...
        'MarkerEdgeColor','black',...
        'MarkerFaceColor','black')
    hold on
end


lim3 = max(defcoorx(:)) + 1;
title('Displacement Field')
xlim([lim2 lim3])
ylim([lim2 lim])


if postchoice == 1
    
    
    %stress field (extra)
    subplot(2,3,4)
    
    for z = 1:NumElem
        for ii = 1:4
            xz(1,ii) = defcoorx(1,C(z,ii));
            yz(1,ii) = defcoory(1,C(z,ii));
            stress1(1,ii) = Se(1,C(z,ii));
        end
        xz(1,5) = defcoorx(1,C(z,1));
        yz(1,5) = defcoory(1,C(z,1));
        plot(xz,yz,"-",'color',[1 0 0],'MarkerSize',3,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        patch(xz,yz,stress1)
        colorbar
        hold on
    end
    
    lim3 = max(defcoorx(:)) + 1;
    title('Stress-x')
    xlim([lim2 lim3])
    ylim([lim2 lim])
    
    %stress field (extra)
    subplot(2,3,5)
    
    for z = 1:NumElem
        for ii = 1:4
            xz(1,ii) = defcoorx(1,C(z,ii));
            yz(1,ii) = defcoory(1,C(z,ii));
            stress1(1,ii) = Se(2,C(z,ii));
        end
        xz(1,5) = defcoorx(1,C(z,1));
        yz(1,5) = defcoory(1,C(z,1));
        plot(xz,yz,"-",'color',[1 0 0],'MarkerSize',3,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        patch(xz,yz,stress1)
        colorbar
        hold on
    end
    
    lim3 = max(defcoorx(:)) + 1;
    title('Stress-y')
    xlim([lim2 lim3])
    ylim([lim2 lim])
    
    %stress field (extra)
    subplot(2,3,6)
    
    for z = 1:NumElem
        for ii = 1:4
            xz(1,ii) = defcoorx(1,C(z,ii));
            yz(1,ii) = defcoory(1,C(z,ii));
            stress1(1,ii) = Se(3,C(z,ii));
        end
        xz(1,5) = defcoorx(1,C(z,1));
        yz(1,5) = defcoory(1,C(z,1));
        plot(xz,yz,"-",'color',[1 0 0],'MarkerSize',3,...
            'MarkerEdgeColor','black',...
            'MarkerFaceColor','black')
        patch(xz,yz,stress1)
        colorbar
        hold on
    end
    
    lim3 = max(defcoorx(:)) + 1;
    title('Stress-xy')
    xlim([lim2 lim3])
    ylim([lim2 lim])
else
end
%%


