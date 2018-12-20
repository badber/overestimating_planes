clear all
close all
clc

%% Inputs
P_Loss_max = 1.8e3;
[P_Deloaded,EFR] = meshgrid(0:10:400,0:10:200);

D = 0.5e-2;
P_D = 45e3;

number_planes = 15;

%% Error check
P_Deload_range = [0 P_Deloaded(end)];
EFR_range = [0 EFR(end)];

if (P_Deload_range(end)-P_Deload_range(1))<(EFR_range(end)-EFR_range(1))
    error('The range of Deload should be greater than the range of EFR, otherwise the plane approximation will not work (I would have to switch Deload and EFR in the calculation of the planes)')
end

%% Conservative fit 
z = @(x,y) (P_Loss_max-x-y).^2*(10/(4*0.8)); % "x" and "y" are P_Deload and EFR, respectively

a=P_Deload_range(end)+EFR_range(end);
segment_diagonal = (a/2*sqrt(2))/number_planes;
clear a
segment = segment_diagonal*sqrt(2);
% To understand how I calculate the value of the segments, just draw the
% domain of Pdeload and EFR and graph how the planes would be calculated.

% Create the points where I will fit planes:

if (P_Deload_range(1)+segment<EFR_range(end))
    EFR_point = P_Deload_range(1)+segment;
    points_plane{1,1} = [P_Deload_range(1), EFR_range(1), z(P_Deload_range(1),EFR_range(1))];
    points_plane{1,2} = [P_Deload_range(1)+segment, EFR_range(1), z(P_Deload_range(1)+segment,EFR_range(1))];
    points_plane{1,3} = [P_Deload_range(1), EFR_point, z(P_Deload_range(1),EFR_point)];
else
    PDeload_point = -1*(EFR_range(end)-(P_Deload_range(1)+segment));
    points_plane{1,1} = [P_Deload_range(1), EFR_range(1), z(P_Deload_range(1),EFR_range(1))];
    points_plane{1,2} = [P_Deload_range(1)+segment, EFR_range(1), z(P_Deload_range(1)+segment,EFR_range(1))];
    points_plane{1,3} = [PDeload_point, EFR_range(end), z(PDeload_point,EFR_range(end))];
end

for i=2:number_planes-1
    if (P_Deload_range(1)+i*segment<EFR_range(end))
        EFR_point = P_Deload_range(1)+i*segment;
        points_plane{i,1} = points_plane{i-1,3};
        points_plane{i,2} = [P_Deload_range(1)+i*segment, EFR_range(1), z(P_Deload_range(1)+i*segment,EFR_range(1))];
        points_plane{i,3} = [P_Deload_range(1), EFR_point, z(P_Deload_range(1),EFR_point)];
    elseif (P_Deload_range(1)+i*segment>P_Deload_range(end))
        EFR_point = P_Deload_range(1)+i*segment-P_Deload_range(end);
        PDeload_point = -1*(EFR_range(end)-(P_Deload_range(1)+i*segment));
        points_plane{i,1} = points_plane{i-1,2};
        points_plane{i,2} = [P_Deload_range(end), EFR_point, z(P_Deload_range(end),EFR_point)];
        points_plane{i,3} = [PDeload_point, EFR_range(end), z(PDeload_point,EFR_range(end))];
    else
        PDeload_point = -1*(EFR_range(end)-(P_Deload_range(1)+i*segment));
        points_plane{i,1} = points_plane{i-1,3};
        points_plane{i,2} = [P_Deload_range(1)+i*segment, EFR_range(1), z(P_Deload_range(1)+i*segment,EFR_range(1))];
        points_plane{i,3} = [PDeload_point, EFR_range(end), z(PDeload_point,EFR_range(end))];
    end
end

% Last plane:
points_plane{number_planes,1} = points_plane{number_planes-1,2};
points_plane{number_planes,2} = points_plane{number_planes-1,3};
points_plane{number_planes,3} = [P_Deload_range(end), EFR_range(end), z(P_Deload_range(end),EFR_range(end))];
   
% Graphically check that the points are properly placed:
fig=1;
figure(fig)
fig=fig+1;
surf(P_Deloaded,EFR,z(P_Deloaded,EFR)*1e-6, 'FaceColor', [255,100,0]/255, 'FaceAlpha', .9);
xlabel('P_{Deload} (MW)')
ylabel('EFR (MW)')
zlabel('$$\frac{\left(\textrm{P}_{\textrm{Loss}}^{\textrm{max}}-P_{\textrm{Deload}}-\textrm{EFR}\right)^2\cdot \textrm{T}_\textrm{d}}{4\cdot \Delta f_{\textrm{max}}} \qquad (\textrm{GW}^2 \cdot \textrm{s}^2$$)','Interpreter','latex')
hold on
view([37.5 20])
for plane=1:size(points_plane,1)
    plot_matrix(1,:) = points_plane{plane,1};
    plot_matrix(2,:) = points_plane{plane,2};
    plot_matrix(3,:) = points_plane{plane,3};
    scatter3(plot_matrix(:,1),plot_matrix(:,2),plot_matrix(:,3)*1e-6,'filled')
end

j=3;
scatter3(plot_matrix(j,1),plot_matrix(j,2),plot_matrix(j,3)*1e-6,'filled')

% Now get the planes defined by those points:
a_vector = [];
b_vector = [];
c_vector = [];
for i=1:number_planes
    vector_plane{i} = cross(points_plane{i,1} - points_plane{i,2},...
        points_plane{i,1} - points_plane{i,3});
    d{i} = points_plane{i,1}*vector_plane{i}';
    
    InputFile_Constants_plane{i} = [vector_plane{i}(1), vector_plane{i}(2), -d{i}]/(-vector_plane{i}(3));
    
    a_vector = [a_vector InputFile_Constants_plane{i}(1)];
    b_vector = [b_vector InputFile_Constants_plane{i}(2)];
    c_vector = [c_vector InputFile_Constants_plane{i}(3)];
    
    % For plotting the planes:
%     Z{i} = ((vector_plane{i}(1)*P_Deloaded +...
%         vector_plane{i}(2)*EFR - d{i})/(-vector_plane{i}(3))...
%         + PFR*EFR*(0.5/(4*0.8)))/PFR;
    
    Z{i} = InputFile_Constants_plane{i}(1)*P_Deloaded +...
        InputFile_Constants_plane{i}(2)*EFR +...
        InputFile_Constants_plane{i}(3);
end

% First plot the theoretical solution:
Z_theoretical = (P_Loss_max-P_Deloaded-EFR).^2*(10/(4*0.8));
figure(fig)
fig=fig+1;
%shading interp
h1 = surf(P_Deloaded,EFR,Z_theoretical*1e-6, 'FaceColor', [255,100,0]/255, 'FaceAlpha', .9);
xlabel('P_{Deload} (MW)')
ylabel('EFR (MW)')
zlabel('$$\frac{\left(\textrm{P}_{\textrm{L}}^{\textrm{max}}-P_{\textrm{DL}}-R_\mathcal{S}\right)^2\cdot \textrm{T}_\textrm{PFR}}{4\cdot \Delta f_{\textrm{max}}} \qquad (\textrm{GW}^2\textrm{s}^2$$)','Interpreter','latex')
%title(['Linearization, ' num2str(number_planes) ' planes'])
hold on
%view([37.5 20])
view([20 25])

zlim([1 11])

% % % For a video:;
% % frame = 1;
% % F(frame) = getframe(gcf);
% % pause(2)
% % frame = frame+1;

% save_fig=1;
% print(figure(fig-1),'-dpng', ['Fig' num2str(save_fig)])
% print(figure(fig-1),'-depsc', ['Fig' num2str(save_fig)])
% save_fig=save_fig+1;

h2 = surf(P_Deloaded,EFR,Z{1}*1e-6, 'FaceColor', [1,255,200]/255, 'FaceAlpha', .7);
legend([h1,h2], {'Exact solution', 'Linear fits'},'Position',[0.65 0.6 0.1 0.14]);

% print(figure(fig-1),'-dpng', ['Fig' num2str(save_fig)])
% print(figure(fig-1),'-depsc', ['Fig' num2str(save_fig)])
% save_fig=save_fig+1;

% % F(frame) = getframe(gcf);
% % pause(0.5)
% % frame = frame+1;
for i=2:number_planes
    surf(P_Deloaded,EFR,Z{i}*1e-6, 'FaceColor', [1,255,200]/255, 'FaceAlpha', .7);
    
%     print(figure(fig-1),'-dpng', ['Fig' num2str(save_fig)])
%     print(figure(fig-1),'-depsc', ['Fig' num2str(save_fig)])
%     save_fig=save_fig+1;
% %         F(frame) = getframe(gcf);
% %         pause(0.5)
% %         frame = frame+1;
end
legend([h1,h2], {'Exact solution', 'Linear fits'},'Position',[0.65 0.6 0.1 0.14]);
title(['Linearization, ' num2str(number_planes) ' planes'])

% %legend([h1,h2], {'Exact solution', 'Linear fits'},'Position',[0.65 0.6 0.1 0.14]);
% % title(['Linearization, ' num2str(number_planes) ' planes'])
% % pause(2)
% % 
% % video = VideoWriter('Marginal.mp4','MPEG-4');
% % video.FrameRate = 1;
% % open(video)
% % writeVideo(video,F)
% % close(video)


%% Check "conservativeness" of the fit:
clear H h1 h2 i PFR plane plot_matrix z Z points_plane vector_plane d

max_P_Deloaded = max(max(max(P_Deloaded,1)));
max_EFR = max(max(max(EFR,1)));
RHS_theoretical = (P_Loss_max-P_Deloaded-EFR).^2*(10/(4*0.8)); % Analytical Right-Hand_side

for i=1:number_planes
    P_Deloaded_plane{i} = zeros(size(P_Deloaded));
    EFR_plane{i} = zeros(size(P_Deloaded));
    RHS_AppropriateArea{i} = zeros(size(P_Deloaded));
    Ones_plane{i} = zeros(size(P_Deloaded));
end
% Get matrices with values only in the points where the plane is approximating
% the curve, not in the rest of the points (put zeros in the rest of the
% points):
for k=1:size(P_Deloaded,1)
    for j=1:size(P_Deloaded,2)
        if (P_Deloaded(k,j)<=-(EFR(k,j)-(P_Deload_range(1)+segment))) 
            P_Deloaded_plane{1}(k,j) = P_Deloaded(k,j);
            EFR_plane{1}(k,j) = EFR(k,j);
            RHS_AppropriateArea{1}(k,j) = RHS_theoretical(k,j);
            Ones_plane{1}(k,j) = 1;
        end
        for index=2:number_planes-1
            if ((P_Deloaded(k,j)>-(EFR(k,j)-(P_Deload_range(1)+segment*(index-1))))...
                    && (P_Deloaded(k,j)<=-(EFR(k,j)-(P_Deload_range(1)+segment*(index)))))
                P_Deloaded_plane{index}(k,j) = P_Deloaded(k,j);
                EFR_plane{index}(k,j) = EFR(k,j);
                RHS_AppropriateArea{index}(k,j) = RHS_theoretical(k,j);
                Ones_plane{index}(k,j) = 1;
            end
        end            
        if (P_Deloaded(k,j)>-(EFR(k,j)-(P_Deload_range(1)+segment*(number_planes-1))))
            P_Deloaded_plane{number_planes}(k,j) = P_Deloaded(k,j);
            EFR_plane{number_planes}(k,j) = EFR(k,j);
            RHS_AppropriateArea{number_planes}(k,j) = RHS_theoretical(k,j);
            Ones_plane{number_planes}(k,j) = 1;
        end
    end
end

for i=1:number_planes
    AppropriateArea_Plane{i} = InputFile_Constants_plane{i}(1)*P_Deloaded_plane{i} + ...
        InputFile_Constants_plane{i}(2)*EFR_plane{i} + ...
        InputFile_Constants_plane{i}(3)*Ones_plane{i};
    differences{i} = AppropriateArea_Plane{i}-RHS_AppropriateArea{i};
    percentages{i} = differences{i}./RHS_AppropriateArea{i}*100;
end

% Maximum overestimation for each plane:
for i=1:number_planes
    max_error_perecentage_plane(i) = max(max(max(percentages{i})));
end
max_error_perecentage_plane
max_ALL_PLANES = max(max_error_perecentage_plane)

% Mean overestimation for all planes:
% (first set to 0 all the NaN entries:)
for i=1:number_planes
    percentages{i}(isnan(percentages{i})) = 0;
end
percentages_total = zeros(size(percentages{1}));
for i=1:number_planes
    percentages_total = percentages_total + percentages{i};
end
mean_error_percentage_ALL_PLANES = mean(mean(percentages_total))

%save('InputFile_Constants_plane')

%% Include the damping term
% First plot the theoretical solution:
Z_theoretical = (P_Loss_max-P_Deloaded-EFR).^2*(10/(4*0.8)) - ((P_Loss_max-P_Deloaded-EFR)*10*D*P_D/4);
figure(fig)
fig=fig+1;
h1 = surf(P_Deloaded,EFR,Z_theoretical*1e-6, 'FaceColor', [255,100,0]/255, 'FaceAlpha', .9);
xlabel('P_{Deload} (MW)')
ylabel('EFR (MW)')
zlabel('H (GW \cdot s^2)')
hold on
view([37.5 20])

% Now get the planes with the damping term
for i=1:number_planes
     
    Z{i} = InputFile_Constants_plane{i}(1)*P_Deloaded +...
        InputFile_Constants_plane{i}(2)*EFR +...
        InputFile_Constants_plane{i}(3)...
        - ((P_Loss_max-P_Deloaded-EFR)*10*D*P_D/4);
end

h2 = surf(P_Deloaded,EFR,Z{1}*1e-6, 'FaceColor', [1,255,200]/255, 'FaceAlpha', .7);
for i=2:number_planes
    surf(P_Deloaded,EFR,Z{i}*1e-6, 'FaceColor', [1,255,200]/255, 'FaceAlpha', .7);
end
legend([h1,h2], {'Exact solution', 'Linear fits'},'Position',[0.65 0.6 0.1 0.14]);
title(['Linearization, ' num2str(number_planes) ' planes'])
