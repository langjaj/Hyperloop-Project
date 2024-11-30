%%Function to calculate the squircle geometry with a fixed blockage ratio
%%given a squareness and aspect ratio.

%% SQUIRCLE EQUATION (x/a)^2 + (y/b)^2 - sqrt(s) (x/a)^2 (y/b)^2 = 1        %%note not s^2 as is often used, this makes it more linear change to square,  use sqrt(s)
%% mainly work with x as a function of y here for ease when exporting to cad
%% equation becomes when using aspect ratio: x = (b/c) * sqrt(  (b^2 - y.^2)  ./  (b^2 - s*y.^2)
%% a is horizontal radius, b is vertical radius, s is squareness, c is aspect ratio b/a

%% Can translate squircle equation down such that the separation distance (lowest point) is always 0.1d=0.3m from the bottom of the pod if the pod is centred at the origin (requires centre to be 1.8m below origin))
%% use change of variables y-->y-b+1.8 i.e. y--> y-ycentre
%% x = (b/c) * sqrt(  (b^2 - (y-b+1.8).^2)  ./  (b^2 - sqrt(s)*(y-b+1.8).^2)


clear all
close all
hold on
delete matlab*coordinates*.txt      %delete all old coordinate files


%%INPUTS
%fixed
beta = 0.25;                      %blockage ratio
nominal_area = (pi*1.5^2)/beta;   %
%variables:
%VECTORS of values  
c_input_calc_below_1 = round([0.4584:(1-0.4584)/4:1],3);  %these values are chosen to roughly evenly space min AR and 1 and then similar between 1/(max AR) and 1, its confusing becasue the AR changes nature above or below 1, this makes sense 
c_input_calc_above_1 = round([1:(1.4057-1)/2:1.4057],3);     %%merge together to get c inputs

c_input = [0.458, 0.594,  0.729,   0.864,  1,  1.203,  1.406];     %[0.459 0.7295 1 1.317];  %0.4584:0.2:1.317981495616843;      %aspect ratio a/b    Cant have additional zeros    
s_input = [0, 0.2, 0.4, 0.6, 0.8, 1];      %squareness


separation_distance = 0.3;   %% this is the distance between the bottom of the pod and the tube
y_low = -1.5-separation_distance;      %%% pod is 1.5m diameter and centred at the origin, so this gives the y coordinate of the bottom of the tube

%write a file of the params to be read into spaceclaim
writematrix(c_input,'aspect_ratio_values.txt')
writematrix(s_input,'squareness_values.txt')




%set colours
newcolors = [219 114 43 ; 90 65 110  ; 188 177 76;  157 188 173 ; 65 144 92;  182 48 36; 52 146 173 ; 152 200 7 ; 201 73 102]/255;
set(0,'DefaultAxesColorOrder',newcolors)

%set all interpreters to latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%set font size
set(groot,'defaultAxesFontSize',8)


%%SQUIRCLE DEFINITION x=....
%squircle_eq_master = @(y,b,c,s) (b./c) * sqrt( (b.^2 - (y).^2)  ./  (b.^2- s.*(y).^2) );   %%%% centred on origin 
squircle_eq_master = @(y,b,c,s) (b./c) * sqrt( (b.^2 - (y-b-y_low).^2)  ./  (b.^2 - (s^(1/2)).*(y-b-y_low).^2) );   %%%% bottom is separation dist from pod from pod, tube centre is (b+y_low)

%%%set up figs
%pod and tubes
fig1 = figure(1);
fig1.Units = 'centimeters';
fig1.Position(3)=9;                               %changes overall figure size in cm [width height]
xlabel('$y$ (m)')
ylabel('$z$ (m)')
xlim([-4.7 4.7])
ylim([-2 5.5])
%set aspect ratio of plot
%daspect([1 1 1]) this messes up sizes so do it manually
ax1 = gca;
ax1.Units = 'centimeters';
ax1.Position(3) = 7.8;     %set width
ax1.Position(4) = ax1.Position(3) /   ( (ax1.XLim(2)-ax1.XLim(1)) / (ax1.YLim(2)-ax1.YLim(1))) ;     %forces correct aspect ratio
ax1.Position(1) = 1.1;         %set right

box on

%design space
fig2 = figure(2);
fig2.Units = 'centimeters';
fig2.Position(3)=9;                               %changes overall figure size in cm [width height]
xlabel('Squareness, $S$')
ylabel('Aspect Ratio, $\mathit{AR} = b/a$')
hold on
xlim([-0.05 1.05])
ylim([0.4 1.5])
ax2 = gca;
ax2.Units = 'centimeters';
ax2.Position = ax1.Position;
ax2.Position = ax1.Position;
box on
xticks(s_input);
yticks(c_input);
xtickformat('%.1f')
ytickformat('%.2f')

%design space but the points are the actual shapes
fig3 = figure(3);
fig3.Units = 'centimeters';
fig3.Position(3)=9;                               %changes overall figure size in cm [width height]
xlabel('Squareness, $S$')
ylabel('Aspect Ratio, $\mathit{AR} = b/a$')
hold on
xlim([-0.05 1.05])
ylim([0.4 1.5])
ax3 = gca;
ax3.Units = 'centimeters';
ax3.Position = ax1.Position;
ax3.Position(4) = ax1.Position(3) /   ( (ax3.XLim(2)-ax3.XLim(1)) / (ax3.YLim(2)-ax3.YLim(1))) ;     %forces correct aspect ratio
box on
xticks(s_input);
yticks(c_input);
xtickformat('%.1f')
ytickformat('%.2f')

%test of max clearance plot
fig4 = figure(4);
hold on


%%%% initialise structure to store x and y values of each shape to export
%%%% and plot later. The array stores and array of [s ar x y] values for
%%%% each (S,AR) combination.
tube_shapes_array = [];

%%CALL FUNCTIONS TO SOLVE
%main cases
v_radius_array = zeros(length(s_input),length(c_input));       %store vertical radius for each case
for i = 1:length(c_input)
    for j = 1:length(s_input)
        v_radius_result = numeric_solve(c_input(i),s_input(j),nominal_area,squircle_eq_master,y_low);
        v_radius_array(j,i) = v_radius_result;
        tube_shapes_array_line = check_and_plot(i,j,c_input(i),s_input(j),v_radius_result,nominal_area,squircle_eq_master,newcolors,y_low);
        tube_shapes_array = [tube_shapes_array  tube_shapes_array_line];     
    end
end
writematrix(tube_shapes_array,'tube_shapes_coordinate_array.txt')    %export array so that it can be plotted in different scripts
%limiting cases
height_limit_case(nominal_area,squircle_eq_master,y_low);
width_limit_case(nominal_area,squircle_eq_master,y_low);


%%POD PLOT
x_pod = -1.5:0.01:1.5;
y_pod = real(sqrt( 1.5^2 - x_pod.^2));
figure(1)
area(x_pod,y_pod,'FaceColor','black',FaceAlpha=0.2);
area(x_pod,-y_pod,'FaceColor','black',FaceAlpha=0.2);


%%%SAVE PLOTS
exportgraphics(figure(1),'pod_and_tubes_for_optimising-matlab.pdf')  
exportgraphics(figure(3),'ar_and_s_tube_design_space.pdf') 
exportgraphics(figure(1),'pod_and_tubes_for_optimising-matlab.png')  
exportgraphics(figure(3),'ar_and_s_tube_design_space.png') 










%%%%%%FUNCTIONS%%%%%%%%%%%%


%DONT USE --- doesnt account for offset etc
function v_radius = symbolic_solve(c,s,y_low)
%SYMBOLIC WAY -- takes way too long for s =/= 0 or 1, dont bother using
%%solve area integral to find horizontal radius
syms y b
integrand = (b/c) * sqrt(  (b^2 - y^2)  /  (b^2 - s*y^2)  );                 %i.e. x = ....
integral = int(integrand,y_pod,[y_low,y_low+b]) ;                                          %integral from bottom (y_low) to centre(y_low+b) (which gives quarter area as only positive brach of equation is being used which is half of geom)
equation = 4 * integral == nominal_area;                                     % *4 to give full area and equate to target area
v_radius_symbolic = abs(vpasolve(equation,b));                                     % solve symbolically
v_radius = vpa(v_radius_symbolic);                                               %convert symbolic answer to number
v_radius = v_radius(1);
end




function v_radius = numeric_solve(c,s,nominal_area,squircle_eq,y_low)
%%NUMERIC WAY
    integrand_def = @(y,b) squircle_eq(y,b,c,s) ;                                        %i.e. x = ....
    integral_def = @(b) integral( @(y) integrand_def(y,b) ,y_low,y_low+b);                          %integral from bottom (y_low) to centre(y_low+b) (which gives quarter area as only positive brach of equation is being used which is half of geom)
    equation = @(b) nominal_area - 4.* integral_def(b);                                   % *4 to give full area and equate to target area
    options = optimoptions('fsolve','OptimalityTolerance',1e-6);
    v_radius = real(fsolve(equation,3.5,options));                                                     %solve equation=0 starting from number (3.5 here)
end






function tube_shapes_array_line = check_and_plot(i,j,c,s,v_radius,nominal_area,squircle_eq,newcolors,y_low)
%%%%check by integrating back again numerically
%%COORDINATES of quarter geometry
if s==1
    y = [-1.8 -1.8 2*v_radius-1.8 2*v_radius-1.8];                        %only needs corner points for a square
    x = [0 v_radius/c v_radius/c 0];
else
    y = -1.8:v_radius/500:2*v_radius-1.8;             %generates 1000 points which is the max desing modeller can take in. All mins are at -1.8 by design 
    x = squircle_eq(y,v_radius,c,s);
end
x=real(x);
%Plot the shapes around the pod:
figure(1)  
plot(x,y,-x,y,'LineWidth',1, 'Color',  [newcolors(i,:) 1])  %(8-j)/8 ])   %change one for this commented bit to get varying transparency     %this goves 3 element vecto for colour and final element is transparancy, the five will need changing if theres more
%Plot the design points:
figure(2)  
scatter(s,c,'x','MarkerEdgeColor',newcolors(i,:),'MarkerEdgeAlpha',(15-j)/15,'LineWidth',2)     %plot design point      ,'MarkerFaceAlpha',0.2
%Plot the desgin points as the actual shapes:
figure(3)     
rf = 150;     %reduction factor to scale each 
fig3xpos = s+(x/rf);
fig3xneg = s+(-x/rf);
fig3y = c+(y-y_low-v_radius)/rf;
plot( fig3xpos,fig3y,fig3xneg,fig3y,'LineWidth',2,'Color', [newcolors(i,:)  (15-j)/15 ])   %change one for this commented bit to get varying transparency
%%%% generate table of x and y values of each shape to export and plot
%%%% later, each array line is is a single cross section, theyre all appended
%%%% together in the loop
tube_shapes_array_line = [s+zeros(1,length(x)) ; c+zeros(1,length(x)) ; x ; y-y_low-v_radius];    %change y to center it at zero again

%%%TEST, calculate maximum distance from the pod (equivalent to maximum
%%%distance from origin - radius), this is just magnitude of (x,y).
figure(4)
dist_from_0 = (x.^2 + y.^2).^(1/2);
max_clearance = max(dist_from_0,[],'all');
scatter(s,max_clearance,'x','MarkerEdgeColor',newcolors(i,:),'MarkerEdgeAlpha',(8-j)/8,'LineWidth',2)


%write to file that can be imported by Spaceclaim
filename = 'matlab_coordinates_AR' + string(c) + '_SQ' + string(s) + '_spaceclaim.txt';
fileID = fopen(filename,'w');          %create and open file with aspect ratio and squareness in filename
%fprintf(fileID,'%12s %12s\n',x,y);
array_format = [-x;y];                                       %make an array with the required format, this stacksx and y vectors, x is negative as this is the right half for the pod geometry
if s==1
    fprintf(fileID,'%s\n%s\n','Polyline=true','3d=true');        %this tells it to create a polyline rather than spline so that it passes the corner poins for s=1
else
    fprintf(fileID,'%s\n%s\n','Polyline=false','3d=true');        %this tells it to create a smooth spline rather than polyline
end
% format: [group(i.e 1)  x-coord y-coord z--coord]
fprintf(fileID,'%12.7f 0 %12.7f\n',array_format);          %  %i is integer, %12.7f is floating point width 12 and 7 after decimal
fprintf(fileID,'\n%12.7f 0 %12.7f',-x(1),y(1),-x(end),y(end));   % add line between first and last points to close shape, this is on new line so its not comibned with the rest which messes up the spline
fclose(fileID);

%check areas with other solution methods
fun = @(y) squircle_eq(y,v_radius,c,s);
area_check2 = 4*integral(fun,y_low,y_low+v_radius) - nominal_area                                %want this to be 0

% % % %write to file that can be imported by DM
% % % filename = 'matlab_coordinates_AR' + string(c) + '_SQ' + string(s) + '_DM.txt';
% % % fileID = fopen(filename,'w');          %create and open file with aspect ratio and squareness in filename
% % % %fprintf(fileID,'%12s %12s\n',x,y);
% % % point_id_counter = 1:length(y);                    %this just numbers each point sequentially
% % % array_format = [point_id_counter;-x;y];               %make an array with the required format 
% % % % format: [group(i.e 1) point_id x-coord y-coord z--coord]
% % % fprintf(fileID,'1 %i 0 %12.7f %12.7f\n',array_format);          %  %i is integer, %12.7f is floating point width 12 and 7 after decimal
% % % fclose(fileID);
end


%SOLVE FOR LIMITING CASE PARAMETERS ENSURING SEPARATION IS ALWAYS > 0.1d =
%0.3m

%HEIGHT LIMIT case for s=1, need to solve for aspect ratio where b = v_rad =
%0.6d = 1.8m
function c_limit_height = height_limit_case(nominal_area,squircle_eq,y_low)
b_limit_1 = 1.8;
s_limit_1= 1;

integrand_def = @(y,c) squircle_eq(y,b_limit_1,c,s_limit_1);
integral_def = @(c) integral( @(y) integrand_def(y,c) ,y_low,y_low+b_limit_1);                          %integral from bottom (y_low) to centre(y_low+b) (which gives quarter area as only positive brach of equation is being used which is half of geom)
integral_equation = @(c) nominal_area - 4.* integral_def(c);                                   % *4 to give full area and equate to target area

options = optimoptions('fsolve','OptimalityTolerance',1e-6);
c_limit_height = fsolve(integral_equation,0.5,options)                                                     %solve equation=0 starting from number (0.5 here)
end



%WIDTH LIMIT case for s=0, need to smiltaneously solve for aspect ratio and height where
%area holds and minimum distance from the origin is dist_limit = 1.8m
function [c_limit_width,b_limit_width] = width_limit_case(nominal_area,squircle_eq,y_low)
s_limit = 0;
dist_limit = 1.8; 

integrand_def = @(y,b,c) squircle_eq(y,b,c,s_limit);                                        %i.e. x = ....
integral_def = @(b,c) integral( @(y) integrand_def(y,b,c) ,y_low,y_low+b);                             %integral from bottom (y_low) to centre(y_low+b) (which gives quarter area as only positive brach of equation is being used which is half of geom)
integral_equation = @(b,c) nominal_area - 4.* integral_def(b,c);                                   % *4 to give full area and equate to target area

distance_fun = @(y,b,c) sqrt( y^2  + squircle_eq(y,b,c,s_limit)^2 );                                    %i.e sqrt(y^2 + x^2), gives distance from the origin
min_dist_y = @(b,c) fminbnd(@(y) distance_fun(y,b,c),-dist_limit,0);                                         %minimisation of th distance, want to find b and c such that this distance is 0, this returns the y value of the minimum
min_dist_eq = @(b,c) distance_fun(min_dist_y(b,c),b,c) - dist_limit;                                                %this evaluates the function at the y value to give the minimum distance (rather than coordinate), this need to be solved to be equalt to 1.8

eq_system = @(X) [integral_equation(X(1),X(2));min_dist_eq(X(1),X(2))];                % package together the two equations, one solving for area and one for minimum distance from the origin
x0 = [1;1];                                                                             %initial point

options = optimoptions('fsolve','OptimalityTolerance',1e-8,'MaxFunEvals',100000);
limits_width = fsolve(eq_system,x0,options);                                            %solve equations simulataneously
b_limit_width = real(limits_width(1))                                                   %split output into b and c and get rid of imaginary part
c_limit_width = real(limits_width(2))
end


%%%% to run spaceclaim from here use
%%%%  !"C:\Program Files\ANSYS Inc\v202\scdm\spaceclaim.exe" /RunScript="D:\OneDrive - University of Leeds\Project_Research\simulations\3d_baseline\optimisation\spaceclaim-generate-script-working.py" /Headless=True /Splash=False /Welcome=False /ExitAfterScript=True






