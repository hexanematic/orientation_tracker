function plot_orientation_field_single_cell(image, CenterOfMass, gamma_vector, p)

% ------ graphical output of orientation field ------
% plot_orientation_field_single_cell(image, CenterOfMass, gamma_vector, p)
%
% Cite this code: Armengol-Collado et al. 2023 (Nature Physics)
%
% Author: Julia Eckert
% Date: 2023-6
%
%
% --- INPUT: 
% image:        2D matrix or image
% CenterOfMass: center of mass of polygons
% gamma_vector: complex shape functions of all polygons
% p:            p-fold order, 
%               2 for nematic or
%               6 for hexatic
%
% example: 
% plot_orientation_field_single_cell(image, CenterOfMass, gamma_vector, 2)
%
%
% --- OUTPUT:
% image ... director field overlayed with topological defects
%           red circles: positive defects
%           blue circles: negative defects
% ------------------------------------------------------

if nargin ~= 4
    disp('Not enough input arguments')
    disp('plot_orientation_field_single_cell(image, CenterOfMass, gamma_vector, p)')
end


% scale factor for length of directors
scale_factor = 10; 

figure
imagesc(image)
axis off
hold on
set(gcf,'color','white')

% center of mass of polygon
xy = CenterOfMass; 

% p-fold shape function
gamma = gamma_vector(:,p);

% calculate angles
if p == 2
    alpha = (angle(gamma))/p;
end
if p == 6
    alpha1 = (angle(gamma))/p;
    alpha2 = (angle(gamma))/p+(pi/(p/2));
    alpha3 = (angle(gamma))/p+(2*pi/(p/2));
    alpha = [alpha1, alpha2, alpha3]; 
end

% plot orientation field
for i = 1:length(alpha(1,:))    
    direc = [cos(alpha(:,i)),sin(alpha(:,i))];

    hDefl = quiver(xy(:,1),xy(:,2), -scale_factor*direc(:,1), -scale_factor*direc(:,2));
    hold on
    set(hDefl,...
        'Color',[1 0.999 1],...
        'LineWidth',2,...
        'MaxHeadSize',0,...
        'AutoScale','off');
    hDefl = quiver(xy(:,1),xy(:,2), scale_factor*direc(:,1), scale_factor*direc(:,2));
    set(hDefl,...
        'Color',[1 0.999 1],...
        'LineWidth',2,...
        'MaxHeadSize',0,...
        'AutoScale','off');
end

set(gcf,'color','white')
set(gcf,'Position',[10 10 500 500])
set(gca,'position',[0 0 1 1],'units','normalized');
set(gca,'Visible','off');
set(gcf,'Units','pixels');
set(gcf,'PaperPositionMode','auto');
saveas(gcf,[num2str(p),'-fold_orientation_field.png'])

end