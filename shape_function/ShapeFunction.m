
function gamma = ShapeFunction(Vertices_all, CenterOfMass)

% ------ p-fold shape function ------
% gamma = ShapeFunction(Vertices_all, CenterOfMass)
%
% Cite this code: Armengol-Collado et al. 2023 (Nature Physics)
%
% Author: Julia Eckert
% Date: 2023-6 
%
%
% --- DESCRIPTION: 
% This code calculates the p-fold shape function, gamma_p, of cells/polygons, 
% as shown in Eq.(2) in the above reference.
%
%
% --- INPUT: 
% This code requires the vertex positions per cell.
%
% Vertices_all: Vertices_all{k,1}(:,1:2) ... structure containing lists of 
%                                            xy vertex positions of cell k
%               ... example: Vertices_all = Polygon.Vertex; 
%
% CenterOfMass: CenterOfMass(:,1:2) ... list of xy positions of cells
%               [] ... xy positions are calculated using vertex positions
%               ... example: CenterOfMass = Polygon.CenterOfMass; 
%
%
% --- OUTPUT:
% gamma.vector(k,p):           complex p-fold shape function of cell k   
%
% gamma.center_of_mass(:,1:2): xy positions of cells
%
% ------------------------------------------------------

if nargin ~= 2
    disp('Not enough input arguments')
    disp('ShapeFunction(Vertices_all,CenterOfMass) or ShapeFunction(Vertices_all,[])')
end


% center of mass of all cell-polygons 
if isempty(CenterOfMass)
    disp('xy positions are calculated using vertex positions)')
    % calculate center of mass per polygon
    Polygon_center_all = [];
    for vl = 1:size(Vertices_all,1)
        Vertex_list = unique(Vertices_all{vl,1}, 'rows');
        x_poly = 1/length(Vertex_list(:,1))*sum(Vertex_list(:,1));
        y_poly = 1/length(Vertex_list(:,2))*sum(Vertex_list(:,2));
        Polygon_center_all = [Polygon_center_all; x_poly, y_poly];
    end
else
    Polygon_center_all = CenterOfMass;
end


% calculate shape function per polygon
gamma_p_vec = [];
for nc = 1:length(Polygon_center_all(:,1))
    
    % center of mass of polygon
    xy_poly = Polygon_center_all(nc,1:2);

    % vertex positions per polygon
    xy_vertex = unique(Vertices_all{nc,1}, 'rows');
       
    % calculate:
    % - distance, r_v, between center of mass (CoM) and vertex 
    % - angle, angles_v, between x-axis and CoM-vertex position
    r_v = [];
    angle_v = [];   
    for nv = 1:length(xy_vertex)
        
        % CoM-vertex distance
        r_nv = abs(norm(xy_vertex(nv,:)' - xy_poly(:)));
        r_v = [r_v; r_nv];
        
        % angle between x-axis and CoM-vertex positions
        angle_nv = atan2(xy_vertex(nv,2)-xy_poly(2), xy_vertex(nv,1)-xy_poly(1));
        angle_v = [angle_v; angle_nv];        
    end
    
    % p-fold shape function
    for p = 1:7
        gamma_p_vec(nc,p) = sum(r_v(:).^p.*exp(1j*p*angle_v(:,1)))/...
                              sum((r_v(:).^p));
    end
end
   
% output:
gamma = [];
gamma.vector = gamma_p_vec;
gamma.center_of_mass = Polygon_center_all;

end



