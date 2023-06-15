### Cite as:
Armengol-Collado et al. 2023 (Nature Physics)

Author: Julia Eckert

Date: 2023-6 


### Scripts:

Each of the following scripts contains a detailed description at the beginning of the script.

```MATLAB
ShapeFunction.m
```
This code computes the p-fold shape function, *gamma_p*, ( for p = 2 ... 6) for each cell in two-dimensional monolayers, as shown in Eq.(2) of the above reference.

```MATLAB
ShapeParameter.m
```
This code coarse-grains the shape functions, *gamma_g*, over a given length scale and computes the shape parameter, *Gamma_p*, as show in Eq.(3) in the above reference.

```MATLAB
plot_orientation_field_single_cell.m
```
This code uses the shape functions, *gamma_p*, to plot the nematic (p=2) and hexatic (p=6) orientation fields of each cell.

```MATLAB
plot_coarse_grained_orientation_field.m
```
This code uses the shape parameters, *Gamma_p*, to plot the coarse-grained nematic (p=2) and hexatic (p=6) orientation field for a given length scale. The orientation field is superimposed with topological defects. 

### EXAMPLE:

```MATLAB
load('Image_Polygon.mat')
```

The Matlab file contains lists of vertex positions **Polygon.Vertex** and center of mass **Polygon.CenterOfMass** per cell *k* :

```MATLAB
Vertices_all = Polygon.Vertex;        %Vertices_all{k,1}(:,1:2) 
CenterOfMass = Polygon.CenterOfMass;  %CenterOfMass(k,1:2)
```

The complex p-fold shape function is obtained by:

```MATLAB
gamma = ShapeFunction(Vertices_all, CenterOfMass)
```

and if the center of mass of the cell is not available:
```MATLAB
gamma = ShapeFunction(Vertices_all, [])
```

Output:
```MATLAB
gamma.vector(k,p)                 % complex p-fold shape function of cell k   
gamma.center_of_mass(:,1:2)       % xy coordinates of cells
```

The complex p-fold shape parameter is determined by coarse-graining over the individual shape functions, **gamma.vector**, within disks of radii, *R*, of integer multiples of *cg_radius* with *R = m * cg_radius*. This code requires the output file **gamma** of **ShapeFunction.m**.

```MATLAB
cg_radius = 20;                   % minimum coarse-graining radius 
grid_distance = 20;               % grid distance to define the coarse-graining disk positions
im = dir(‘Image.tif');            
data = imread(im.name);
image_size = size(data);          % image size for coarse-graining
```
```MATLAB
Gamma = ShapeParameter(gamma, cg_radius, grid_distance, image_size)
```
Output:
```MATLAB
Gamma.xy_disk(:,1:2)      % xy coordinates of disks in grid configuration of distance cg_radius
Gamma.R                   % radii of disks = coarse-graining radii
Gamma.vector{l}(R,p)      % complex p-fold shape parameter in each grid point l for all radii R
Gamma.angle{l}(R,p)       % p-fold orientation of shape parameter in each grid point l for all radii R
```
