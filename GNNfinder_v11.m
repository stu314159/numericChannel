%% Version 11
clear
close all
clc
% This version requires the geom2d toolbox to be installed.
% If flip = 1, ensure that the airfoil points are initially defined in a counterclockwise direction. If
% they are not, make sq_frac negative.

% add path to geom2d toolbox
addpath('./geom2d','./geom2d/geom2d');
addpath('./geom2d/polygons2d');
addpath('./geom2d/polynomialCurves2d');


%% CHOOSE GEOMETRY
fprintf('loading the airfoil geometry data.\n');
% Airfoil
A = importdata('eee_blade5_scaled.dat');  % airfoil data
% Interpolate the outline
A = interppolygon(A,100);
A = interppolygon(A,100); % twice to ensure even spacing

%% INPUTS
% Flip Airfoil?
flip = 1; 
    % Flip the airfoil upside down if flip = 1
    if flip == 1;
       A(:,2) = -A(:,2);
    else
    end
% Lattice Inputs
Ny_divs = 65;
Lx_p = 4;
Ly_p = 3;
Lz_p = 14;
AoA = 45; % stagger angle of attack (measured from -z axis up to the chord line), degrees
ceiling = 0.85*Ly_p;
% Squealer Inputs
sq_on = 1;      % 1 to turn on squealer
sq_frac = 1/20; % Squealer fraction - fraction of Lo that represents the distance the squealer is pushed in from the outside tip edge
    % it should never be less than 1/25
sq_depth = 1/6; % Squealer depth - fraction of the blade depth that the squealer removes

%% MAXIMIZE DISTANCE BETWEEN POINTS ON THE BLADE
distance = 0;
point1 = [];
point2 = [];
for i = 1:length(A)
  x1(i) = A(i,1);
  y1(i) = A(i,2);
  for j = 1:length(A)
      x2(j) = A(j,1);
      y2(j) = A(j,2);
      distance_to_test = sqrt((x1(i)-x2(j))^2 + (y1(i)-y2(j))^2);
      if distance_to_test > distance
          distance = distance_to_test;
          i_point1 = i;
          j_point2 = j;
          point1 = [x1(i),y1(i)];
          point2 = [x2(j),y2(j)];
      else
      end
  end
end

%% FIND THE LINE THAT CONTAINS THE TWO POINTS OF MAX DISTANCE
% Create the max distance line
m_dist = (point1(2) - point2(2))/(point1(1) - point2(1));
x_dist = linspace(-10,10,100);
b_dist = point1(2) - m_dist*point1(1);
y_dist = m_dist*x_dist + b_dist;
% Find its angle to the x axis
theta_dist = atan(m_dist); 

%% ROTATE BLADE SO THE POINTS OF MAX DISTANCE ALIGN WITH THE X AXIS
% Keep old x and y values for calculation
A_old = A;
% Perform Rotation
A(:,1) = A_old(:,1).*cos(-theta_dist) - A_old(:,2).*sin(-theta_dist);
A(:,2) = A_old(:,1).*sin(-theta_dist) + A_old(:,2).*cos(-theta_dist);
% Find point 1 and 2 again
point1 = [A(i_point1,1),A(i_point1,2)]; % note that the indices of A have not changed after rotation
point2 = [A(j_point2,1),A(j_point2,2)];

%% TRANSLATE BLADE SO THE POINTS OF MAX DISTANCE ALIGN WITH THE X AXIS
% Find distance to translate vertically
d2t_v = -point1(2);
% Perform vertical translation
A(:,2) = A(:,2) + d2t_v;
% Find distance to translate horizontally 
    %(so that all points have positive x value)
d2t_h = min(A(:,1));
if d2t_h > 0
    d2t_h = 0; % don't translate if all already have positive x values
end
% Perform horizontal translation
A(:,1) = A(:,1) + d2t_h;    
% Find point 1 and 2 again
point1 = [A(i_point1,1),A(i_point1,2)]; % note that the indices of A have not changed after rotation, so this will work
point2 = [A(j_point2,1),A(j_point2,2)];

%% FIND HIGHEST AND LOWEST POINTS FROM X AXIS
% Maximize, minimize y values
    % This works because the max distance line is aligned with the x axis
[highest_y,highest_index] = max(A(:,2));
highest_point = A(highest_index,:);
[lowest_y,lowest_index] = min(A(:,2));
lowest_point = A(lowest_index,:);

%% CREATE CHORD LINE
% Find the points for the beginning and end of the camber line
% Trailing edge point (the one furthest from the y axis)
fprintf('creating chord line.\n');
point1dist = point1(1);
point2dist = point2(1);
if point2dist > point1dist
    trail_point = point2;
    lead_point = point1;
else
    trail_point = point1;
    lead_point = point2;
end
% Find index of the trail point and lead point
i_trail_point = find(A(:,1) == trail_point(1) & A(:,2) == trail_point(2));
i_lead_point = find(A(:,1) == lead_point(1) & A(:,2) == lead_point(2));
% Point 3: whichever high/low point has the greater distance from the
% trailing edge
high_point_dist = sqrt((highest_point(:,1) - trail_point(1))^2 + (highest_point(:,2) - trail_point(2))^2);
low_point_dist = sqrt((lowest_point(:,1) - trail_point(1))^2 + (lowest_point(:,2) - trail_point(2))^2);
if high_point_dist > low_point_dist
    point3 = highest_point;
else
    point3 = lowest_point;
end
% Find index of point 3
i_point3 = find(A(:,1) == point3(1) & A(:,2) == point3(2));
% Define the chord line, plot the chord line
m_chord = (point3(2) - trail_point(2))/(point3(1) - trail_point(1));
x_chord = linspace(-10,10,100);
b_chord = point3(2) - m_chord*point3(1);
y_chord = m_chord*x_chord + b_chord;
% Save old A values
A_old = A;
% Save old chord line values
x_chord_old = x_chord;
y_chord_old = y_chord;
% Save old trail point and point 3 values
point3_old = point3;
trail_point_old = trail_point;

%% FIND THE CHARACTERISTIC LENGTH: CHORD LINE 
Lo = sqrt((point3(1)-trail_point(2))^2 + (point3(2)-trail_point(2))^2);

%% CONSTRUCT GCOORD XZ GRID
fprintf('constructing gcoord array \n');
xm = 0;
ym = 0;
zm = 0;
xp = Lx_p;
yp = Ly_p;
zp = Lz_p;
Ny = ceil((Ny_divs - 1)*(Ly_p/Lo))+1;
Nx = ceil((Ny_divs - 1)*(Lx_p/Lo))+1;
Nz = ceil((Ny_divs - 1)*(Lz_p/Lo))+1;
x_space = linspace(xm,xp,Nx);
y_space = linspace(ym,yp,Ny);
z_space = linspace(zm,zp,Nz);
%[X,Y,Z]=meshgrid(x_space,y_space,z_space);
gcoord = zeros((Nx*Ny*Nz),3);
nd = 0;
for z = 1:Nz
    for y = 1:Ny
        for x = 1:Nx
            nd = nd+1;
            gcoord(nd,1)=x_space(x);
            gcoord(nd,2)=y_space(y);
            gcoord(nd,3)=z_space(z);
        end
    end
end

%% SCALE THE BLADE TO DESIRED SIZE
A(:,1) = A(:,1)/(2*Lo);
A(:,2) = A(:,2)/(2*Lo);

%% TRANSLATE AIRFOIL CoM TO ORIGIN
% Find airfoil CoM
z_bar = (sum(A(:,1)))/(length(A(:,1)));
x_bar = (sum(A(:,2)))/(length(A(:,2)));
% Move x dimension down by x_bar
A(:,2) = A(:,2) - x_bar;
% Move z dimension over by z_bar
A(:,1) = A(:,1) - z_bar;

%% RE-DETERMINE AIRFOIL CoM AFTER MOVEMENT
% Airfoil
z_bar = (sum(A(:,1)))/(length(A(:,1)));
x_bar = (sum(A(:,2)))/(length(A(:,2)));

%% ROTATE AND TRANSLATE THE AIRFOIL
% Rotate to desired AoA
A_old2 = A;
current_AoA = atan(m_chord);
desired_AoA = -AoA*pi/180;
AoA = (desired_AoA - current_AoA); % radians
A(:,1) = A_old2(:,1).*cos(AoA) - A_old2(:,2).*sin(AoA);
A(:,2) = A_old2(:,1).*sin(AoA) + A_old2(:,2).*cos(AoA);
% Translate to center of channel
A(:,1) = A(:,1) + 0.4*Lz_p;
A(:,2) = A(:,2) + 0.6*Lx_p;

%% DEFINE THE SQUEALER
fprintf('defining the squealer.\n');
if sq_on == 1
    % Interpolate the outline
    A = interppolygon(A,100);
    A = interppolygon(A,100); % twice to ensure even spacing
    % Cut out the squealer
    B = expandPolygon(A,sq_frac*Lo);
    B = B{1};
    % Remove NaNs, infs
    B(any(isnan(B),2),:)=[]; % remove NaNs
    B(any(isinf(B),2),:)=[]; % remove infs
    % Find important intersection points
        % Find points where A intersects B
        [zi,xi] = polyxpoly(A(:,1),A(:,2),B(:,1),B(:,2));
        % Points where B intersects itself
            % Break B into 2 segments
            B1 = B((1:floor(length(B)/2)),:);
            B2 = B((ceil(length(B)/2)):(length(B)),:);
            [zi2,xi2] = polyxpoly(B1(:,1),B1(:,2),B2(:,1),B2(:,2),'unique');
        % Remove point that is shared in both B arrays
        samepoint = intersect(B1,B2,'rows');
        i = find(zi2==samepoint(:,1) & xi2==samepoint(:,2));
        zi2(i) = [];
        xi2(i) = [];
    % Combined matrix of all intersection points
    allpoints = [[zi,xi];[zi2,xi2]];
    % Only perform next steps if there are 3 intersections
    if length(allpoints) == 3 
        % Determine bounds of points to remove
        zmin = min(allpoints(:,1));
        zmax = max(allpoints(:,1));
        xmin = min(allpoints(:,2));
        xmax = max(allpoints(:,2));
        % For loop to identify and remove the points in the overlap area
        i = 1;
        overlap_points = [];
        indices = [];
        for i = 1:length(B)
            zval = B(i,1);
            xval = B(i,2);
            if (zval>zmin) && (zval<zmax) && (xval>xmin) && (xval<xmax)
               overlap_points = [overlap_points;[zval,xval]];
               indices = [indices;i];
            else
            end
        end
        % Remove the overlap points
        zcol_B = B(:,1);
        xcol_B = B(:,2);
        zcol_B(indices) = [];
        xcol_B(indices) = [];
        B = [zcol_B,xcol_B];
        % Find indices of the points of B that are outside/on A
        i=1;
        outside_matrix = [];
        indices = [];
        for i = 1:length(B)
            [in,on] = inpolygon(B(i,1),B(i,2),A(:,1),A(:,2));
            if in == 0 || on == 1
               outside_matrix = [outside_matrix;[B(i,1),B(i,2)]];
               indices = [indices;i];
            else
            end
        end
        % Remove the outside points
        zcol_B = B(:,1);
        xcol_B = B(:,2);
        zcol_B(indices) = [];
        xcol_B(indices) = [];
        B = [zcol_B,xcol_B];
    % Only perform next steps if there are < 3 intersections
    elseif (length(allpoints)< 3) & (length(allpoints)~=0)
        % find point of max dist from xi2, yi2
        j = 1;
        distance = 0;
        point3_sq = [];
        for j = 1:length(B)
          z3_sq(j) = B(j,1);
          x3_sq(j) = B(j,2);
          distance_to_test = sqrt((zi2(i)-z3_sq(j))^2 + (xi2(i)-x3_sq(j))^2);
          if distance_to_test > distance
              distance = distance_to_test;
              j_point3_sq = j;
              point3_sq = [z3_sq(j),x3_sq(j)];
          else
          end
        end
        % Define rectangle for the polygon containing all inside squealer points
        twopoints = [[zi2,xi2];point3_sq];
        zmin = min(twopoints(:,1));
        zmax = max(twopoints(:,1));
        xmin = min(twopoints(:,2));
        xmax = max(twopoints(:,2));
        rectangle = [[zmin,xmin];[zmin,xmax];[zmax,xmax];[zmax,xmin];[zmin,xmin]];
        rectangle = interppolygon(rectangle,100);
        % Find indices of outside points
        i=1;
        outside_matrix = [];
        indices = [];
        for i = 1:length(B)
            [in,on] = inpolygon(B(i,1),B(i,2),rectangle(:,1),rectangle(:,2));
            if in == 0
               outside_matrix = [outside_matrix;[B(i,1),B(i,2)]];
               indices = [indices;i];
            else
            end
        end
        % Remove the outside points
        zcol_B = B(:,1);
        xcol_B = B(:,2);
        zcol_B(indices) = [];
        xcol_B(indices) = [];
        B = [zcol_B,xcol_B];
    else
    end
  % Interpolate B to get evenly spaced points
  B = interppolygon(B,50);
else
end

%% DETERMINE INSIDE AND OUTSIDE POINTS FOR THE AIRFOIL
% Create the test grid from gcoord
% z values to test
fprintf('Determine inside and outside points for the airfoil\n')
zzlim1 = find((gcoord(:,3))>=(min(A(:,1))));
zzlim2 = find((gcoord(:,3))<=(max(A(:,1))));
zzlims = intersect(zzlim1,zzlim2);
zz = gcoord(zzlims,3);
z_gcoord = unique(zz); % remove duplicate z values
% x values to test
xxlim1 = find((gcoord(:,1))>=(min(A(:,2))));
xxlim2 = find((gcoord(:,1))<=(max(A(:,2))));
xxlims = intersect(xxlim1,xxlim2);
xx = gcoord(xxlims,1);
x_gcoord = unique(xx); % remove duplicate x values
% Test if inside A
i = 1;
j = 1;
inside_matrix = [];
for i = 1:length(z_gcoord)
    for j = 1:length(x_gcoord)
        [in,on] = inpolygon(z_gcoord(i),x_gcoord(j),A(:,1),A(:,2));
        if in == 1 || on == 1
            inside_matrix = [inside_matrix;[z_gcoord(i),x_gcoord(j)]];
        else
        end
        in = 0;
    end
end

%% DETERMINE INSIDE AND OUTSIDE POINTS FOR THE SQUEALER
fprintf('determine inside and outside points for the squealer.\n');
if sq_on == 1
    % Keep the same test grid as the airfoil test points (x_gcoord,z_gcoord)
    % Test if inside B
    i = 1;
    j = 1;
    inside_matrix_sq = [];
    for i = 1:length(z_gcoord)
        for j = 1:length(x_gcoord)
            [in,on] = inpolygon(z_gcoord(i),x_gcoord(j),B(:,1),B(:,2));
            if in == 1 || on == 1
                inside_matrix_sq = [inside_matrix_sq;[z_gcoord(i),x_gcoord(j)]];
            else
            end
            in = 0;
        end
    end
else
end
%% DETERMINE LIST OF AIRFOIL GLOBAL NODE NUMBERS TO ADD TO OBSTACLE FILE
% Change inside_matrix to match gcoord dimensions (x,y,z columns)
fprintf('getting global node numbers for the airfoil.\n');
Y = zeros(length(inside_matrix),1);
inside_matrix = [inside_matrix(:,2),Y,inside_matrix(:,1)]; % Airfoil
% Create test grid (already have x and z values for the interior points)
% y values to test
yylim1 = find((gcoord(:,2))>0);
yylim2 = find((gcoord(:,2))<=ceiling);
yylims = intersect(yylim1,yylim2);
yy = gcoord(yylims,2);
y_gcoord = unique(yy); % remove duplicate y values
% Create new inside_matrix including all y planes
inside_matrix_add = [];
for i=1:length(y_gcoord)
    for j=1:length(inside_matrix)
        inside_matrix_add = [inside_matrix_add;[inside_matrix(j,1),y_gcoord(i,1),inside_matrix(j,3)]];
    end
end
inside_matrix = [inside_matrix;inside_matrix_add];

%% %% DETERMINE LIST OF SQUEALER GLOBAL NODE NUMBERS TO DELETE FROM OBSTACLE FILE
fprintf('getting global node numbers for the squealer.\n');
if sq_on == 1
    % Change inside_matrix to match gcoord dimensions (x,y,z columns)
    Y_sq = zeros(length(inside_matrix_sq),1);
    inside_matrix_sq = [inside_matrix_sq(:,2),Y_sq,inside_matrix_sq(:,1)]; % Squealer
    % Create test grid (already have x and z values for the interior points)
    % y values to test
    yylim1_sq = find((gcoord(:,2))>=(ceiling-sq_depth*ceiling));
    yylim2_sq = find((gcoord(:,2))<=ceiling);
    yylims_sq = intersect(yylim1_sq,yylim2_sq);
    yy_sq = gcoord(yylims_sq,2);
    y_gcoord_sq = unique(yy_sq); % remove duplicate y values
    % Create new inside_matrix including all y planes
    inside_matrix_add_sq = [];
    for i=1:length(y_gcoord_sq)
        for j=1:length(inside_matrix_sq)
            inside_matrix_add_sq = [inside_matrix_add_sq;[inside_matrix_sq(j,1),y_gcoord_sq(i,1),inside_matrix_sq(j,3)]];
        end
    end
    inside_matrix_sq = []; % get rid of the points on the y=0 plane
    inside_matrix_sq = [inside_matrix_sq;inside_matrix_add_sq];
else
end
%% Find the gnn for each inside_matrix (airfoil) lattice point
fprintf('Finding global node numbers for each interior lattice point.\n')
% tic;
% gnn = [];
% i=1;
% for i = 1:length(inside_matrix)
%     newgnn = find(gcoord(:,1)==inside_matrix(i,1) & gcoord(:,2)==inside_matrix(i,2) & gcoord(:,3)==inside_matrix(i,3));
%     gnn = [gnn;newgnn];
% end
% if sq_on == 1
%     % Find the gnn for each inside_matrix_sq (squealer) lattice point
%     gnn_sq = [];
%     i=1;
%     for i = 1:length(inside_matrix_sq)
%         newgnn_sq = find(gcoord(:,1)==inside_matrix_sq(i,1) & gcoord(:,2)==inside_matrix_sq(i,2) & gcoord(:,3)==inside_matrix_sq(i,3));
%         gnn_sq = [gnn_sq;newgnn_sq];
%     end
% else
% end
% time_check = toc;
% 
% fprintf('time for first method = %g.\n',time_check);


% gnn = sort(gnn);
% gnn_sq = sort(gnn_sq);
% 
%tic;
gnn = find(ismember(gcoord,inside_matrix,'rows'));
gnn_sq = find(ismember(gcoord,inside_matrix_sq,'rows'));
%time_check = toc;

%fprintf('time for the second method = %g.\n',time_check);



%fprintf('Plotting the results.\n');
% Plot the solid nodes with global node numbers in the gnn array
%figure(1)
% % Plot blade nodes in black
% i=1;
% for i=1:length(gnn)
%     scatter3(gcoord((gnn(i,1)),1),gcoord((gnn(i,1)),2),gcoord((gnn(i,1)),3),'filled','k')
%     hold on
% end
% hold on
% % Plot squealer nodes in red
% if sq_on == 1
%     i=1;
%     for i=1:length(gnn_sq)
%         scatter3(gcoord((gnn_sq(i,1)),1),gcoord((gnn_sq(i,1)),2),gcoord((gnn_sq(i,1)),3),'filled','r')
%         hold on
%     end
%     hold off
% else
% end

%scatter3(gcoord(gnn,1),gcoord(gnn,2),gcoord(gnn,3),'filled','k');
%if sq_on == 1
%    hold on
%    scatter3(gcoord(gnn_sq,1),gcoord(gnn_sq,2),gcoord(gnn_sq,3),'filled','r');
%    hold off
%end


if sq_on == 1
    % Remove the squealer gnns
    gnn = setxor(gnn,gnn_sq);
else
end

%% save gnn to file

turb_blade_data = 'turbine_blade.mat';
save(turb_blade_data,'gnn','Ny_divs','Lo','Lx_p','Ly_p','Lz_p');
