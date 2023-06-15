
function Gamma = ShapeParameter(gamma, cg_radius, grid_distance, image_size)

% ------ p-fold shape parameter ------
% Gamma = ShapeParameter(gamma, cg_radius, grid_distance, image_size)
%
% Cite this code: Armengol-Collado et al. 2023 (Nature Physics)
%
% Author: Julia Eckert
% Date: 2023-6
%
%
% --- DESCRIPTION: 
% This code coarse-grains the shape functions, gamma_p, (see Eq.(2) in 
% the above reference) over a length scale, R, of integer multiples of 
% cg_radius with R = m * cg_radius. 
% The resulting p-fold shape parameter, Gamma_p, is shown in Eq.(3) 
% in the above reference.
%
%
% --- INPUT: 
% This code requiers the output file of ShapeFunction.mat
%
% gamma:        output file of ShapeFunction.mat; 
%
% cg_radius:    minimum coarse-graning radius and bin-size for n*cg_radius 
%               ... example: cg_radius = 20;
%
% grid_distance: grid size of coarse-graning grid  
%               ... example: grid_distance = 20;
%
% image_size:   size of image 
%               ... example: image_size = [512,512] = size(image);
%
%
% --- OUTPUT:
% Gamma.xy_disk(:,1:2): xy coordinates of disks in grid configuration of distance cg_radius
% Gamma.R:              radii of disks = coarse-graining radii
% Gamma.vector{l}(R,p): complex p-fold shape parameter in each grid point l for all radii R
% Gamma.angle{l}(R,p):  p-fold orientation of shape parameter in each grid point l for all radii R
% ------------------------------------------------------

if nargin ~= 4
    disp('Not enough input arguments')
    disp('ShapeParameter(gamma, cg_radius, grid_distance, image_size)')
end


% input variables 
r_grid = grid_distance;
Polygon_center_all = gamma.center_of_mass;
gamma_p_vec = gamma.vector;
 
% grid with distance r_grid
xi = r_grid:r_grid:image_size(2)-r_grid;
yi = r_grid:r_grid:image_size(1)-r_grid;

grid_points = [];
for i = 1:length(xi)
    for j = 1:length(yi)
        grid_points = [grid_points; xi(i), yi(j)];
    end
end

% disk size binning
center_point_n = [image_size(1),image_size(2)]/2;
bin = [cg_radius:cg_radius:center_point_n(:,1)];

% calculate shape parameter for all grid points
for g = 1:length(grid_points(:,1))    
    % distance between grid point and center of mass of polygons
    dis_grid_poly = [];
    for k = 1:length(Polygon_center_all(:,1))
        xy_poly = Polygon_center_all(k,1:2);
        dis_grid_poly = [dis_grid_poly; abs(norm(xy_poly-grid_points(g,:)))];  
    end

    % calculate p-fold shape parameter for each disk size
    glob_mean_gamma_vec = [];
    glob_mean_gamma_angle = [];
    for p = 1:7
        for d = 1:length(bin)
            dis = bin(d);
            % complex gamma shape functions in disk of radius bin
            gamma_bin = gamma_p_vec(dis_grid_poly<=dis,p);
                        
            % mean of complex shape functions
            mean_gamma_bin = nanmean(gamma_bin);
            glob_mean_gamma_vec(d,p) = mean_gamma_bin;

            % angle of shape parameter
            glob_mean_gamma_angle(d,p) = angle(mean_gamma_bin);

        end
    end
    % p-fold shape parameter per grid-point g
    Gamma_p_vec{g,1}(:,:) = glob_mean_gamma_vec;    
    % orientation of the p-fold shape parameter per grid-point g
    Gamma_p_angle{g,1}(:,:) = glob_mean_gamma_angle;
end

% output:
Gamma = [];
Gamma.angle = Gamma_p_angle;
Gamma.vector = Gamma_p_vec;
Gamma.xy_disk = grid_points;
Gamma.R = bin;
        
end



