function defect_cg = plot_coarse_grained_orientation_field(image, Gamma, p, cg)

% ------ coarse-grained orientation field + topological defects ------
% defect_cg = plot_coarse_grained_orientation_field(image, Gamma, p, cg)
%
% Cite this code: Armengol-Collado et al. 2023 (Nature Physics)
%
% Author: Julia Eckert
% Date: 2023-5
%
%
% --- INPUT: 
% image:        2D matrix or image
% Gamma:        output file of ShapeParameter.mat; 
% p:            p-fold order, 
%               2 for nematic or
%               6 for hexatic
% cg:           coarse graining step-size of Gamma.R (1,2,...N)
%
% example: 
% plot_coarse_grained_orientation_field(image, Gamma, 2, 4)
%
%
% --- OUTPUT:
% defect_cg.charge:  topological charge
% defect_cg.xy:      xy position of defects
% defect_cg.p:       p-fold order
% defect_cg.cg:      coarse-graining step (1,2,...N)
% image ... director field overlayed with topological defects
%           red circles: positive defects
%           blue circles: negatice defects
% ------------------------------------------------------

if nargin ~= 4
    disp('Not enough input arguments')
    disp('plot_coarse_grained_orientation_field(image, Gamma, p, cg)')
end


% scale factor for length of directors
scale_factor = 10; 

% input variables 
Gamma_p_vec = Gamma.vector;
xy = Gamma.xy_disk;

figure
imagesc(image)
hold on

beta = [];
for h = 1:size(Gamma_p_vec,1) 

    if p == 2
        beta = [beta; angle(Gamma_p_vec{h,1}(cg,p))/p];
    end
    if p == 6
        beta1 = angle(Gamma_p_vec{h,1}(cg,p))/p;
        beta2 = beta1+(pi/(p/2));
        beta3 = beta1+(2*pi/(p/2));
        beta = [beta; beta1, beta2, beta3]; 
    end

end

for i = 1:length(beta(1,:))

    direc = [cos(beta(:,i)),sin(beta(:,i))];

    hDefl = quiver(xy(:,1), xy(:,2), -scale_factor*direc(:,1), -scale_factor*direc(:,2));
    hold on
    set(hDefl,...
        'Color',[1 0.999 1],...
        'LineWidth',2,...
        'MaxHeadSize',0,...
        'AutoScale','off');
    hDefl=quiver(xy(:,1), xy(:,2), scale_factor*direc(:,1), scale_factor*direc(:,2));
    hold on
    set(hDefl,...
        'Color',[1 0.999 1],...
        'LineWidth',2,...
        'MaxHeadSize',0,...
        'AutoScale','off');
end

%%.............defect charge using shape parameter...........
[uni_val_x,uni_int_x] = unique(xy(:,1),'first');
[uni_val_y,uni_int_y] = unique(xy(:,2),'first');

defect_charge = [];
xy_defect = [];

% get squared director configuration and calculate defect charge
int_beta = beta;
int_x = xy(:,1);
int_y = xy(:,2);

for uiy = 1:length(uni_val_y)-1
    uiy = uiy-1;
    for uix = 1:length(uni_val_x)-1   
        list4 = [];
        xylist4 = [];
        if uni_int_x(uix)+1 <= length(int_y)
            list4 = [int_beta(uni_int_x(uix)+uiy), int_beta(uni_int_x(uix+1)+uiy),...
                int_beta(uni_int_x(uix+1)+1+uiy), int_beta(uni_int_x(uix)+1+uiy),...
                int_beta(uni_int_x(uix)+uiy)];  
            xylist4 = [int_x(uni_int_x(uix)+uiy), int_y(uni_int_x(uix)+uiy); ...
                int_x(uni_int_x(uix+1)+uiy), int_y(uni_int_x(uix+1)+uiy); ...
                int_x(uni_int_x(uix+1)+1+uiy), int_y(uni_int_x(uix+1)+1+uiy); ...
                int_x(uni_int_x(uix)+1+uiy), int_y(uni_int_x(uix)+1+uiy)];
        end
        % position of defect
        polyin = polyshape(xylist4(:,1),xylist4(:,2));
        [xwind,ywind] = centroid(polyin);

        angle2_diff = [];
        for j = 1:length(list4)-1
            % angle difference of director orientation
            angle2_dif = list4(j+1)-list4(j);
            % winding number
            if angle2_dif > pi/p
                angle2_dif = angle2_dif - 2*pi/p;
            end
            if angle2_dif < -pi/p
                angle2_dif = angle2_dif + 2*pi/p;
            end
            angle2_diff = [angle2_diff; angle2_dif];
        end
        % defect charge   
        defect_charge = [defect_charge; sum(angle2_diff)/(2*pi)];
        % position defect charge
        xy_defect = [xy_defect; xwind, ywind];
    end
end


hold on
XY = xy_defect;
DC = defect_charge;
if p == 2
    for i = 1:length(XY(:,1))
    
        if DC(i) < -0.8 && DC(i) > -1.2 %-1
            plot(XY(i,1), XY(i,2),'o', 'color',[0 0.4470 0.7410], 'MarkerSize', 10, 'LineWidth', 2)
        end
        if DC(i) < -0.3 && DC(i) > -0.7 %-1/2
            plot(XY(i,1), XY(i,2), 'o','color',[0.3010 0.7450 0.9330], 'MarkerSize', 10, 'LineWidth', 2)
        end
        if DC(i) < 0.7 && DC(i) > 0.3 %+1/2
            plot(XY(i,1), XY(i,2),'o', 'color',[0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 2)
        end    
        if DC(i) < 1.2 && DC(i) > 0.8 %+1
            plot(XY(i,1), XY(i,2),'o', 'color',[0.6350 0.0780 0.1840], 'MarkerSize', 10, 'LineWidth', 2)
        end
    end
end
if p == 6
    for i=1:length(XY(:,1))
        
        if DC(i) < -0.3 && DC(i) > -0.4 %-2/6
            plot(XY(i,1), XY(i,2),'o', 'color',[0 0.4470 0.7410], 'MarkerSize', 10, 'LineWidth', 2)
        end
        if DC(i) < -0.1 && DC(i) > -0.2 %-1/6
            plot(XY(i,1), XY(i,2), 'o','color',[0.3010 0.7450 0.9330], 'MarkerSize', 10, 'LineWidth', 2)
        end
        if DC(i) < 0.2 && DC(i) > 0.1 %+1/6
            plot(XY(i,1), XY(i,2),'o', 'color',[0.8500 0.3250 0.0980], 'MarkerSize', 10, 'LineWidth', 2)
        end    
        if DC(i) < 0.4 && DC(i) > 0.3 %+2/6
            plot(XY(i,1), XY(i,2),'o', 'color',[0.6350 0.0780 0.1840], 'MarkerSize', 10, 'LineWidth', 2)
        end
    end
end

set(gcf,'color','white')
set(gcf,'Position',[10 10 500 500])
set(gca,'position',[0 0 1 1],'units','normalized');
set(gca,'Visible','off');
set(gcf,'Units','pixels');
set(gcf,'PaperPositionMode','auto');
saveas(gcf,[num2str(p),'-fold_orientation_field_cg-step_',num2str(cg),'.png'])

defect_cg = [];
defect_cg.charge = defect_charge;
defect_cg.xy = xy_defect;
defect_cg.p = p;
defect_cg.cg = cg;

end