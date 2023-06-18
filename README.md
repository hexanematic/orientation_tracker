### Cite as:
Armengol-Collado *et al.*, "Epithelia are multiscale active liquid crystals", Nat. Phys. (2023).

### Data:
The folder 'data' contains the data support the findings of this study. 
The subfolders are organized as follows.

Folder 'dat' contains text file named 'exp_#.dat', with # a number ranging from 1 to 68 (corresponding to each analyzed experimental configuration).
Each line in these files refers to a cell in the correspondig microscopy picture.
The first column contains the x-coordinate of the centre of mass of the cell. 
The second column contains the y-coordinate of the centre of mass of the cell. 
The third column contains the real part of the complex shape function gamma_2 of the cell. 
The fourth column contains the imaginary part of the complex shape function gamma_2 of the cell. 
The fifth column contains the real part of the complex shape function gamma_3 of the cell and so on.

Folder 'mat' contain Matlab files 'Image#_Polygon.mat' contains the centre of mass per cell and the vertices per cell, for the corresponding image.

Folder 'tif' contains experimental microscopy pictures named 'image#.tif'

Folder 'xlsx' contains spreadsheets with the source data presented in Armengol-Collado et al. 2023 (Nature Physics).

### Scripts:

Each of the following scripts contains a detailed description at the beginning of the script.

```MATLAB
ShapeFunction.m
```
This script computes the *p*-fold shape function of epithelial monolayers, defined as
```math
\gamma_p = \frac{\sum_v |\boldsymbol{r}_{v}|^{p}e^{ip\phi_{v}}}{\sum_v |\boldsymbol{r}_v|^p}
```
where $`\boldsymbol{r}_{v}=\{x_{v},y_{v}\}`$ and $`\phi_{v}=\arctan(y_{v}/x_{v})`$, with $`v=1,\,2\ldots\,V`$, are the position and orientation of the $`v-`$th vertex of a $`V-`$sided polygon representing the contour of a segmented cell. 

```MATLAB
ShapeParameter.m
```
This script compute the *p*-fold shape parameter of epithelial monolayers, defined upon coarse-graining the shape function $`\gamma_p`$ over the length scale $`R`$. That is
```math
\Gamma_p = \langle \gamma_p \rangle_R 
```
where the average $`\langle\cdots\rangle_R`$ denotes an average over all the cells whose center lies within a distance $`R`$ from the point at which $`\Gamma_p`$ is computed.

```MATLAB
plot_orientation_field_single_cell.m
```
This code uses the shape functions, $`\gamma_p`$, to plot the nematic ($`p=2`$) and hexatic ($`p=6`$) orientation fields of each cell.

```MATLAB
plot_coarse_grained_orientation_field.m
```
This code uses the shape parameters, $`\Gamma_p`$, to plot the coarse-grained nematic ($`p=2`$) and hexatic ($`p=6`$) orientation field for a given length scale. The orientation field is superimposed with topological defects. 

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
