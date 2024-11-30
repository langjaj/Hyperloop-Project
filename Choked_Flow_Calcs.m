clear all
close all
clc

%blockage ratio
beta = 0.25;

%specific heat ratio and specific gas constant
g=1.4;
R = 287.05;

%limit equations
syms M
isentropic_eqn = beta == 1 -( M .* ( ((g+1)./2) .^ ( (g+1)./(2*(g-1)) ) )  .*   (   (1+ (((g-1)/2) .* (M.^2))  ).^(-(g+1)/(2*(g-1))   )   )   );
kantrowitz_eqn = beta == 1 - M.* ((( (g+1) * (M.^2)) ./ ((g - 1) * (M.^2) + 2)).^(-g/(g - 1))) .* (((g + 1)./(2 * g * (M.^2) - (g - 1))).^(-1/(g - 1))) .* ((1 + ((g - 1)/2)* (M.^2)).^(-(g + 1)./(2*(g - 1)))) .* (((g + 1)/2).^((g + 1)/(2*(g - 1)))) ;

%solve for limits
isentropic_limit_sub = vpasolve(isentropic_eqn,M,[0,1])
isentropic_limit_super = vpasolve(isentropic_eqn,M,[1,50])
kantrowitz_limit = vpasolve(kantrowitz_eqn,M)

Msup = isentropic_limit_super;
attach_eqn2 = M^2  ==  (1 + ((g - 1)/2)* (Msup.^2))  ./  ( g*(Msup.^2)  - ((g - 1)/2) ); 
attach_limit2 = vpasolve(attach_eqn2,M,[0,2])


global beta g R isentropic_limit_sub isentropic_limit_super kantrowitz_limit


%ISENTROPIC FLOW 
 
%INPUTS!!!
%tube temp
T_tube = 300;
%tube pressure
p_tube = 100;
rho_tube = p_tube/(R*T_tube);

%single pod Mach number (relative to tube temperature, absolute coords),
%outputs this to console%
M_pod = 1.1;
[T_isen_sub,T_isen_super,p_isen_sub,u_isen_sub,p_isen_super,choke_strength,u_isen_super,v_lead_shock_podcoords,M_shock_down,v_shock_down_podcoords] = calc(M_pod,T_tube,p_tube)



%%%array of Mach Numbers -- outputs to graphs%%%
M_pod_array = 0.5:0.1:3;
T_isen_sub_array = zeros(1,length(M_pod_array));
T_isen_super_array = zeros(1,length(M_pod_array));
p_isen_sub_array = zeros(1,length(M_pod_array));
u_isen_array = zeros(1,length(M_pod_array));
p_isen_super_array = zeros(1,length(M_pod_array));
choke_strength_array = zeros(1,length(M_pod_array));
for i = 1:length(M_pod_array)
    [T_isen_sub_array(i),T_isen_super_array(i),p_isen_sub_array(i),u_isen_array(i),p_isen_super_array(i),choke_strength_array(i)] = calc(M_pod_array(i),T_tube,p_tube);
end

% %these were a test of whether something proprtional to the drag coefficient
% %could be prdicted by computing the force due to the upstream and
% %downstream pressure difference (i.e. equivalent to asuuming the pod is a rectangle
%v_pod_array = M_pod_array .* sqrt(g*R*T_tube);
%dragish_array = (p_isen_sub_array-p_isen_super_array)*3.14*1.5^2% ./ (0.5*rho_tube*(v_pod_array.^2));

%%%choked simulation results
M_pod_array_sims = [ 0.55 0.6 0.7 0.8 1.1 1.7 2 2.5 2.8];
T_isen_sub_array_sims = [305.77148 311.23105 322.38681 333.74011  368.9816 444.78296 485.34216 557.3269 602.61322];
p_isen_array_sims = [6.9103775 13.754485 28.704386 45.310066 106.61089 297.79468 440.1864 777.27319 1053.7046] +100 ;   %add 100 as values in fluent are gauge pressure
u_isen_array_sims = [174.27759 176.02126 179.25247 182.48039  192.157 211.47911 221.12904 236.32118 245.99475];
M_isen_array_sims = [0.49734709 0.49789774 0.49818772 0.49845785  0.49919519  0.50039148  0.50088531 0.49953297 0.50006151];


%%%PLOTTING%%%
newcolors = {'#db722b','#9dbcad','#bcb14c','#5a416e','#b63024','#41905c','#3492ad','#98c807','#c94966'};
newcolors = [219 114 43 ; 157 188 173 ; 188 177 76 ; 90 65 110 ; 182 48 36 ; 65 144 92; 52 146 173 ; 152 200 7; 201 73 102]/255;
set(0,'DefaultAxesColorOrder',newcolors)
set(0,'DefaultLineLineWidth',1)
set(0,'DefaultAxesFontSize',8)
%set all interpreters to latex
list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

fig = figure('units','centimeters');
fig.Position(3)=9;                               %changes overall figure size in cm [width height]
fig.Position(4) = 12;
colororder(newcolors)

ax = axes(fig,'visible','off');                             %set up common x/y label, needs to go after subplot set up 
ax.XLabel.Visible='on';
ax.YLabel.Visible='on';
ylab = ylabel('Upstream Choked Values');
%ylab.Position(1) = -0.115;

subplot_y_position_variable = 0.35;

%velocity plot
ax1 = axes('Position',[0.13 0.73 0.78 0.2]);
axis_set_up()
plot(M_pod_array,u_isen_array,'Color',newcolors(1,:))
plot(M_pod_array_sims,u_isen_array_sims,'x','MarkerEdgeColor','black')
ylab = ylabel('Velocity (m/s)');
ylab.Position(1) = subplot_y_position_variable;
ylim([145 310])
xticklabels({})

%pressure plot
ax2 = axes('Position',[0.13 0.52 0.78 0.2]);
axis_set_up()
plot(M_pod_array,p_isen_sub_array,'Color',newcolors(2,:))
plot(M_pod_array_sims,p_isen_array_sims,'x','MarkerEdgeColor','black')
ylab = ylabel('Pressure (Pa)');
ylab.Position(1) = subplot_y_position_variable;
ylab.Position(2) = 720;
ylim([0 1400])
yticks([0 400 800 1200])
xticklabels({})

%temperature plot
ax3 = axes('Position',[0.13 0.31 0.78 0.2]);
axis_set_up()
plot(M_pod_array,T_isen_sub_array,'Color',newcolors(3,:))
plot(M_pod_array_sims,T_isen_sub_array_sims,'x','MarkerEdgeColor','black')
ylab = ylabel('Temperature (K)');
ylab.Position(1) = subplot_y_position_variable;
ylab.Position(2) = 520;
ylim([170 890])
xticklabels({})

%machno plot
ax4 = axes('Position',[0.13 0.1 0.78 0.2]);
axis_set_up()
isentropic_line = yline(0.5034,'-','Color',newcolors(4,:),'LineWidth',1,'DisplayName','1D Isentropic Theory');
plot(M_pod_array_sims,M_isen_array_sims,'x','MarkerEdgeColor','black','DisplayName','3D Simulation')
xlabel('Pod Mach Number')
ylab = ylabel('Mach Number');
ylab.Position(1) = subplot_y_position_variable;
ylab.Position(2) = 0.493;
%ylab.Position(2) = -50;
ylim([0.45 0.53])
yticks([0.46 0.48 0.5 0.52])
xticklabels('auto')
legend('Location','southeast')

%add margin to stop croppping too tight. can be adjusted with numbers [bottomleftX , bottomleftY, width, height]
a = annotation('rectangle',[0.03 0.03 0.955 0.91],'Color','w');
exportgraphics(gcf,'choked-flow-theory-simulation-comparison-beta0.25-axes.pdf') 
delete(a)





%CHOKED STRENGTH FIG

fig = figure('units','centimeters');
fig.Position(3)=9;                               %changes overall figure size in cm [width height]
fig.Position(4) = 12;
colororder(newcolors)

axis_set_up()
plot(M_pod_array,choke_strength_array,'Color',newcolors(3,:))
%plot(M_pod_array_sims,M_isen_array_sims,'x','MarkerEdgeColor','black','DisplayName','3D Simulation')
xlabel('Pod Mach Number')
ylab = ylabel('Choking Strength');
ylab.Position(1) = subplot_y_position_variable;
ylab.Position(2) = 0.493;
%ylab.Position(2) = -50;
%ylim([0.45 0.53])
%yticks([0.46 0.48 0.5 0.52])
xticklabels('auto')
legend('Location','southeast')










%%%axis set up function 
function y = axis_set_up()
hold on
grid off
box on
title('')
%xlim([1000000 30000000])
ax = gca;
ax.YAxisLocation = 'right';
%set(gca, 'XScale', 'log')
%xticks([1000000 2000000 5000000 10000000 15000000 23000000])
end



%%%calculation_function%%

function [T_isen_sub,T_isen_super,p_isen_sub,u_isen_sub,p_isen_super,choke_strength,u_isen_super,v_lead_shock_podcoords,M_shock_down,v_shock_down_podcoords] = calc(M_pod,T_tube,p_tube)

    global beta g R isentropic_limit_sub isentropic_limit_super kantrowitz_limit

    %%%UPSTREAM SHOCK SPEED AND ISENTROPIC CONDITIONS

    %pod speed (absolute coords), this is equal to upstream speed in pod
    %coords
    v_pod = M_pod*sqrt(g*R*T_tube);
    
    %upstream shock Mach number (absolute coords, relative to tube temp)
    syms msh M
    lead_shock_equation = isentropic_limit_sub == (msh*( (g+1)*M_pod - 2*msh) +2)  /  sqrt( (2*g*(msh^2) -g+1) * ( (g-1)*(msh^2) +2)  );
    M_lead_shock = vpasolve(lead_shock_equation,msh,[0.8,]);       %2 branches of solution, need the positive sonic one
    
    %upstream shock speed (absolute coords, equivalent to upstream speed in shock coords)
    v_lead_shock = M_lead_shock*sqrt(g*R*T_tube);
    
    %upstream shock speed (pod coords, separation speed)
    v_lead_shock_podcoords = v_lead_shock - v_pod;
    
    %isentropic values behind normal shock (calculated in shock coords)
    T_isen_sub = T_tube * (1+ (2*g/(g+1))*((M_lead_shock^2)-1) ) * (2+ ((g-1)*(M_lead_shock^2)) ) / ( (g+1)*M_lead_shock^2 );
    p_isen_sub = p_tube * ((2*g*(M_lead_shock^2)) -g+1)  /  (g+1) ;
    p_isen_sub = p_tube * (1+ (2*g/(g+1)) *(-1+M_lead_shock^2) ) ;
    rho_isen_sub = p_isen_sub/(R*T_isen_sub);
    u_isen_sub_shockcoords = v_lead_shock * ( ( (2+ (g-1)*M_lead_shock^2) ) / ( ((g+1)*(M_lead_shock^2)) ) );


    %isentropic speed in pod fixed coords
    u_isen_sub = u_isen_sub_shockcoords + v_pod - v_lead_shock;
    

    %SHOCK ATTACH DETACH CALCULATION - supersonic temp comes from
    %isentropic relations of temp as flow is assumed to accelerate
    %isentopically around pod -- therefore do Tsuper over Tsub equation
    T_isen_super = T_isen_sub * (2+(g-1)*isentropic_limit_sub^2)   /   (2+(g-1)*isentropic_limit_super^2);
    p_isen_super = p_isen_sub * ( (T_isen_sub/T_isen_super) ^ (-g/(g-1)) );
    u_isen_super = isentropic_limit_super*sqrt(g*R*T_isen_super);

    syms m_isen_s msh_down 
    M_isentropic_shock_fixed_eq = m_isen_s == isentropic_limit_super + sqrt(T_tube/T_isen_super) * (msh_down - M_pod);
    downstream_shock_equation = msh_down^2 == ( 2+(g-1)*m_isen_s^2 )  /  (2*g*(m_isen_s^2) -g+1) ;
    %sub in top equation to bottom rather than trying to solve together as
    %this doesnt='t seem to work
    downstream_shock_equation2 = msh_down^2 == ( 2+(g-1)*(isentropic_limit_super + sqrt(T_tube/T_isen_super) * (msh_down - M_pod))^2 )  /  (2*g*((isentropic_limit_super + sqrt(T_tube/T_isen_super) * (msh_down - M_pod))^2) -g+1);
    
    
    %%%M_shock_down = vpasolve([M_isentropic_shock_fixed_eq,downstream_shock_equation],msh_down,[0 Inf]);
    M_shock_down = vpasolve(downstream_shock_equation2 ,msh_down,[0 Inf]);
    
    %downstream shock speed (absolute coords, equivalent to downstream speed in shock coords)
    v_shock_down = M_shock_down*sqrt(g*R*T_tube);
    
    %downstream shock speed (pod coords, separation speed)
    v_shock_down_podcoords =   v_shock_down - v_pod;


    %%%CHOKING STRENGTH CALC --- i.e mass flow rate upstream vs mass flow
    %%%rate through throat
    
    %calculate conditions upstream of shock in pod coords
    v_upstream = v_pod;
    rho_upstream = p_tube/(R*T_tube);

    %calculate conditions at throat
    T_throat = T_isen_sub * (2+ (g-1)*(isentropic_limit_sub^2)) / (1+g);      %comes rom applying eq 8.40  anderson (2.3.6 in proj) at throat and behind shock
    u_throat = sqrt(g*R*T_throat);
    rho_throat = rho_isen_sub * ( (T_throat/T_isen_sub)^(1/(g-1)) );

    %choke strength -- i.e. ratio of upstream mass flow to throat mass flow
    choke_strength = ((1-beta)) * (rho_throat*u_throat)/(rho_upstream*v_upstream) ;
    choke_strength = double(1/choke_strength);   %convert to number

end
