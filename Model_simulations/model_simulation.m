%% Model simulation

%MODEL_SIMULATION implements models (1) and (4) using a finite element
%method. Note that the PDEToolbox is required.

% Author: Lukas Eigentler (leigentler001@dundee.ac.uk)
% License: GNU GPL
% Last updated: 21/06/2021



clear; 
close all;

%% parameters
global d kmax kcap c Cxloc_in_B1 Cyloc_in_B1 Cxloc_in_B2 Cyloc_in_B2 N time_start sigma
d = 1; kmax = 1; kcap = 10; c = 0.2; % model parameters; set kcap = 1 and c = 0 for isogenic strains
sigma = 1; % switch between isogenic (sigma = 0) and non isogenic (sigma = 1) strains
r = [0,10*sqrt(50)]; % radius of spatial domain

%% initialise pde
N = 2; %no of equations
model = createpde(N);

%% geometry
NC_in = 6; % number of initial microcolonies
NC_in_B1 = 0.5; % ratio of those allocated as B1
NC_in_B1 = ceil(NC_in_B1*NC_in);
NC_in_B2 = NC_in-NC_in_B1;

R1 = [1;0;0;r(end)]; 
sf = 'R001';
ns = char('R001');
gd = R1;
ns = ns';
[dl,bt] = decsg(gd,sf,ns);
geometryFromEdges(model,dl);

%% boundary conditions
applyBoundaryCondition(model,'neumann','edge',1:4,'g',zeros(N),'q',zeros(N));

%% PDE coefficients
specifyCoefficients(model,'f',@fcoeffunction_nondim,'m',0,'d',1,'c',@ccoeffunction_nondim,'a',0);

%% mesh
mesh = generateMesh(model, 'hmax',0.125*sqrt(50),'hmin',0.125*sqrt(50), 'GeometricOrder', 'Linear');

%% initial condition
point_range = 2*sqrt(50); %radius of \Omega_0 in which initial microcolonies are placed

for kk = 1:length(model.Mesh.Nodes(1,:))
    nrm(kk) = norm(model.Mesh.Nodes(:,kk));
end
ic_ind = find(nrm<point_range);
ic_nodes =  model.Mesh.Nodes(:,ic_ind);
for ff = 1:NC_in
    select_ind = randi(length(ic_nodes(1,:)));
    entry = ic_nodes(:,select_ind);
    Cxloc_in(ff) = entry(1);
    Cyloc_in(ff) = entry(2);
    ic_nodes(:,select_ind) = [];

end
Cxloc_in_B1 = Cxloc_in(1:NC_in_B1); Cyloc_in_B1 = Cyloc_in(1:NC_in_B1);
Cxloc_in_B2 = Cxloc_in(NC_in_B1+1:end); Cyloc_in_B2 = Cyloc_in(NC_in_B1+1:end);


% u0both_non_uniform = @u0fun_both_non_uniform; % spatially homogeneous IC
% setInitialConditions(model,u0eq);

u0both_non_uniform = @u0fun_both_non_uniform; % microcolonies IC
setInitialConditions(model,u0both_non_uniform); 

%% solve
time_start = tic;
tlist = linspace(0,25,10); % vector with output times
try
    result = solvepde(model,tlist); % for time-dependent problems
catch
    disp('Solver failure!')
end
toc(time_start)

B10 = result.NodalSolution(:,1,1);
B20 = result.NodalSolution(:,2,1);
B1 = result.NodalSolution(:,1,end);
B2 = result.NodalSolution(:,2,end);

%% Check if solution is feasible
x = linspace(-r(end),r(end),50);
[X,Y] = meshgrid(x,x);
interpol_sol = interpolateSolution(result,X,Y,[1,2],length(result.SolutionTimes));
interpol_sol_B1 = reshape(interpol_sol(:,1),size(X));
interpol_sol_B2 = reshape(interpol_sol(:,2),size(X));
interpol_sol_B2_neg = interpol_sol_B2; interpol_sol_B2_neg(interpol_sol_B2>0) = NaN;
interpol_sol_B1_neg = interpol_sol_B1; interpol_sol_B1_neg(interpol_sol_B1>0) = NaN;
biofilm_ind = find(interpol_sol_B1+interpol_sol_B2>0.5*kcap);
if min(interpol_sol_B1(biofilm_ind))<-2e-2*kcap || min(interpol_sol_B2(biofilm_ind)) <-2e-2*kcap
disp('Warning: infeasible solution detected')
end

%% Calculate access to free space score
B1_ind = find(B10 ~=0); % find B1 locations
B2_ind = find(B20 ~=0); % find B2 locations
B_ind = [B1_ind;B2_ind];
tm = length(B_ind); %total mass
com = sum(result.Mesh.Nodes(:,B_ind),2)/tm; %centre of mass
dt = 2*pi/1000; % Step size in iteration through circle
theta_vec = dt:dt:2*pi; % initialise circle
loop_length = length(theta_vec);
R_in = 10; % radius of circle
initial_theta_b1 = 0; % initialise access to free space score
initial_theta_b1_loc = zeros(1,loop_length); % initialise log of B1 circle segments (for plotting only)
for tt_loop = 1:loop_length % iterate through circle
    tt = theta_vec(tt_loop);
    x = com(1) + R_in*cos(tt); % convert to cartesian coord.
    y = com(2) + R_in*sin(tt);
    cart_coord = [x;y];
    dist_B1 = zeros(1,length(B1_ind));
    dist_B2 = zeros(1,length(B2_ind));                
    for cc = 1:length(B1_ind) % find distances to all B1 microcolonies
        dist_B1(cc) = norm(result.Mesh.Nodes(:,B1_ind(cc)) - cart_coord); 
    end
    for cc = 1:length(B2_ind) % find distances to all B2 microcolonies
        dist_B2(cc) = norm(result.Mesh.Nodes(:,B2_ind(cc)) - cart_coord);
    end  
    [dist_B1_min,dist_B1_min_ind]  = min(dist_B1); % find distance to closest B1 microcolony
    [dist_B2_min,dist_B2_min_ind]  = min(dist_B2); % find distance to closest B2 microcolony
    [~,min_ind] = min([dist_B1_min,dist_B2_min]); % determine which one is closer

    if min_ind == 1 % if B1 is closest, then add to access to free space score
        initial_theta_b1_loc(tt_loop) = 1;
        initial_theta_b1 = initial_theta_b1 +   dt;                  
    end
end
initial_theta_b1 = initial_theta_b1/(2*pi); % normalise to unity

%% plot outcome
f = figure(10);
mapsize = 100;
map1 = [0:1/mapsize:1; zeros(1,mapsize+1); 0:1/mapsize:1]';
map2 = [zeros(1,mapsize+1); 0:1/mapsize:1; zeros(1,mapsize+1)]';
Axes(3) = axes('Position', [0.8/3+0.1 0.05 0.8/3 0.8]); 
pdeplot(model,'XYData',B1,'Contour','off','ColorBar','off','FaceAlpha',1)
hold on
pdeplot(model,'XYData',0*B2,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
xlabel('$B_1$','interpreter','latex')
pbaspect([1 1 1])
caxis manual
caxis([0 kcap]);

Axes(4) = axes('Position', [2*0.8/3+0.15 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',0*B1,'Contour','off','FaceAlpha',1,'ColorBar','off')
hold on
pdeplot(model,'XYData',B2,'Contour','off','ColorBar','off','FaceAlpha',0.5)
xlabel('$B_2$','interpreter','latex')
pbaspect([1 1 1])
caxis manual
caxis([0 kcap]);

Axes(1) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pbaspect([1 1 1])
hold on
pdeplot(model,'XYData',B1,'Contour','off','FaceAlpha',1,'ColorBar','off')
xlabel('Merged','interpreter','latex')
caxis manual
caxis([0 kcap]);

Axes(2) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',B2,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
hold on
hold on
pbaspect([1 1 1])
caxis manual
caxis([0 kcap]);
plot([-2.5; 2.5], [-65,-65], '-w', 'LineWidth', 2)

set(Axes(2), 'XTick', []);
set(Axes(2), 'YTick', []);
set(Axes(2), 'visible', 'off');
set(Axes(1), 'XTick', []);
set(Axes(1), 'YTick', []);
set(Axes(3), 'XTick', []);
set(Axes(3), 'YTick', []);
set(Axes(4), 'XTick', []);
set(Axes(4), 'YTick', []);

colormap(Axes(1),map1)
colormap(Axes(2),map2)
colormap(Axes(3),map1)
colormap(Axes(4),map2)

set(f,'Windowstyle','normal')
set(findall(f,'-property','FontSize'),'FontSize',11)
set(f,'Units','centimeters')
set(f,'Position',[10 10 3*17/4 5])

%% plot initial condition
clear Axes
R_in = 3*sqrt(50);
f1 = figure(11);
Axes(3) = axes('Position', [0.8/3+0.1 0.05 0.8/3 0.8]); 
pdeplot(model,'XYData',B10,'Contour','off','ColorBar','off','FaceAlpha',1)
hold on
pdeplot(model,'XYData',0*B20,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
xlabel('$B_1$','interpreter','latex')
pbaspect([1 1 1])
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);

Axes(4) = axes('Position', [2*0.8/3+0.15 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',0*B10,'Contour','off','FaceAlpha',1,'ColorBar','off')
hold on
pdeplot(model,'XYData',B20,'Contour','off','ColorBar','off','FaceAlpha',0.5)
xlabel('$B_2$','interpreter','latex')
pbaspect([1 1 1])
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);

Axes(1) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pbaspect([1 1 1])
hold on
pdeplot(model,'XYData',B10,'Contour','off','FaceAlpha',1,'ColorBar','off')
xlabel('Merged','interpreter','latex')
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
caxis manual
caxis([0 kcap]);

Axes(2) = axes('Position', [0.05 0.05 0.8/3 0.8]);
pdeplot(model,'XYData',B20,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
hold on
plot(com(1) + R_in*cos(theta_vec(initial_theta_b1_loc(:)==1)),com(2) + R_in*sin(theta_vec(initial_theta_b1_loc(:)==1)), '.', 'color', [1, 0, 1])
plot(com(1) + R_in*cos(theta_vec(initial_theta_b1_loc(:)==0)),com(2) + R_in*sin(theta_vec(initial_theta_b1_loc(:)==0)), '.', 'color', [0, 1, 0])
pbaspect([1 1 1])
ylim([com(2)-R_in,com(2)+R_in])
xlim([com(1)-R_in,com(1)+R_in])
title(['$AFS_1 = ', num2str(initial_theta_b1), '$'], 'interpreter','latex')
caxis manual
caxis([0 kcap]);
plot([com(1)+R_in-6; com(1)+R_in-1], [com(2)-R_in+1; com(2)-R_in+1], '-w', 'LineWidth', 2)

set(Axes(2), 'XTick', []);
set(Axes(2), 'YTick', []);
set(Axes(2), 'visible', 'off');
set(Axes(1), 'XTick', []);
set(Axes(1), 'YTick', []);
set(Axes(3), 'XTick', []);
set(Axes(3), 'YTick', []);
set(Axes(4), 'XTick', []);
set(Axes(4), 'YTick', []);

colormap(Axes(1),map1)
colormap(Axes(2),map2)
colormap(Axes(3),map1)
colormap(Axes(4),map2)

set(f1,'Windowstyle','normal')
set(findall(f1,'-property','FontSize'),'FontSize',11)
set(f1,'Units','centimeters')
set(f1,'Position',[10 10 3*17/4 5])
%% video
% video_ask = input('Create a video of the solution? (Y/n)', 's');
% if video_ask == 'Y' || video_ask == 'y'
%     video_name = input('Enter video filename','s');
%     video_name_print = strrep(video_name,'_',' ');
%     figure
% v = VideoWriter([video_name,'.avi']);
% v.FrameRate = 10;
% open(v);
% 
% for plotind = 1:length(tlist) 
% B1 = result.NodalSolution(:,1,plotind);
% B2 = result.NodalSolution(:,2,plotind);
% close all
% Axes(3) = axes('Position', [0.8/3+0.1 0.05 0.8/3 0.8]); 
% pdeplot(model,'XYData',B1,'Contour','off','ColorBar','off','FaceAlpha',1)
% hold on
% pdeplot(model,'XYData',0*B2,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
% sgtitle('(b) low founder density', 'interpreter','latex')
% xlabel('$B_1$','interpreter','latex')
% % caxis manual
% % caxis([0 kcap]);
% pbaspect([1 1 1])
% 
% Axes(4) = axes('Position', [2*0.8/3+0.15 0.05 0.8/3 0.8]);
% pdeplot(model,'XYData',0*B1,'Contour','off','FaceAlpha',1,'ColorBar','off')
% hold on
% pdeplot(model,'XYData',B2,'Contour','off','ColorBar','off','FaceAlpha',0.5)
% xlabel('$B_2$','interpreter','latex')
% % caxis manual
% % caxis([0 kcap]);
% pbaspect([1 1 1])
% 
% Axes(1) = axes('Position', [0.05 0.05 0.8/3 0.8]);
% pbaspect([1 1 1])
% hold on
% pdeplot(model,'XYData',B1,'Contour','off','FaceAlpha',1,'ColorBar','off')
% xlabel('Merged','interpreter','latex')
% % caxis manual
% % caxis([0 kcap]);
% Axes(2) = axes('Position', [0.05 0.05 0.8/3 0.8]);
% 
% 
% pdeplot(model,'XYData',B2,'Contour','off','FaceAlpha',0.5,'ColorBar','off')
% % caxis manual
% % caxis([0 kcap]);
% pbaspect([1 1 1])
% 
% 
% set(Axes(2), 'XTick', []);
% set(Axes(2), 'YTick', []);
% set(Axes(2), 'visible', 'off');
% set(Axes(1), 'XTick', []);
% set(Axes(1), 'YTick', []);
% set(Axes(3), 'XTick', []);
% set(Axes(3), 'YTick', []);
% set(Axes(4), 'XTick', []);
% set(Axes(4), 'YTick', []);
% 
% colormap(Axes(1),map1)
% colormap(Axes(2),map2)
% colormap(Axes(3),map1)
% colormap(Axes(4),map2)
% 
% 
% sgtitle({['$t = ',num2str(tlist(plotind)), ', k_{max}=',num2str(kmax),', k_{cap}=',num2str(kcap),', c_{12}=',num2str(c12),', c_{21}=',num2str(c21),', d_{B_1}=',num2str(db1), ', d_{B_2}=',num2str(db2),'$'],video_name_print}, 'interpreter','latex')
% 
% frame = getframe(gcf);
% writeVideo(v,frame);
% end
% 
% close(v);
% end

