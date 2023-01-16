function [A,small_diagonals,large_diagonals] = GridAdjacencyMatrix3D(Nx,Ny,Nz)
%GRIDADJACENCYMATRIX adjacency matrix of all 8 neighbors in 2D-grid.
%
% NOTE THE CODE IS INSPIRED BY 
%
% https://stackoverflow.com/questions/3277541/construct-adjacency-matrix-in-matlab/3283732#3283732
% 
% Create the adjacency matrix with all neighbors connected.
%
% Example of such grid:
%    / X / X  +  
%   +---+---+/|
%   | X | X | +
%   +---+---+/|
%   | X | X | +
%   +---+---+/ Hard to draw but multiple connected 
% Plusses are particles lines, dashes and X's are connections.
%
% INPUT
%
%   Nx - (int) Number of particles in x direction.
%
%   Ny - (int) Number of particles in y direction.
%
%   Nz - (int) Number of particles in z direction.
%
% OUTPUT
%
%   A - (mat) The adjacency matrix of all connected points.
%
%   diagonals - (mat) The points of diagonal springs.
%
NP = Nx*Ny*Nz;
x_direction = repmat([ones(Ny*(Nx-1),1);zeros(Ny,1)],Nz,1);
x_direction = x_direction(1:end-Ny);
y_direction = repmat([ones(Ny-1, 1); 0], Nx*Nz, 1);
% Remove last
y_direction = y_direction(1:end-1);
z_direction = repmat(ones(Ny*Nx, 1), Nz-1, 1);
% z_direction = z_direction(1:end-1);
% Diagonals withinplane xy
bl_to_tr_xy = repmat([ones(Ny-1, 1); 0], Nx-1, 1);
bl_to_tr_xy = repmat([bl_to_tr_xy;zeros(Ny,1)],Nz,1);
bl_to_tr_xy = bl_to_tr_xy(1:end-Ny-1);

tl_to_br_xy = repmat([0;ones(Ny-1, 1)], Nx-1, 1);
tl_to_br_xy = repmat([tl_to_br_xy;zeros(Ny,1)],Nz,1);
tl_to_br_xy = tl_to_br_xy(1:end-Ny+1);
% Diagonals xz
tl_to_br_xz = repmat([ones(Ny*(Nx-1),1);zeros(Ny,1)],Nz-1,1);
tl_to_br_xz = tl_to_br_xz(1:end-Ny);
bl_to_tr_xz = repmat([zeros(Ny,1);ones(Ny*(Nx-1),1)],Nz-1,1);
bl_to_tr_xz = [bl_to_tr_xz;zeros(Ny,1)];
% Diagonals yz
tr_to_bl_yz = repmat([ones(Ny-1,1); 0],Nx*(Nz-1),1);
tr_to_bl_yz = tr_to_bl_yz(1:end-1);

tl_to_br_yz = repmat([0;ones(Ny-1,1)],Nx*(Nz-1),1);
tl_to_br_yz(end+1) = 0;

% multidimensional diagonals
% Numbers indicate a cube with 1-4 clockwise on the top
% and 5-8 clockwise on the bottom.
% Example figure:
%    
%     2----3
%    / |  /|
%   1----4 |
%   |  6-|-7
%   | /  |/
%   5----8
%
one_to_seven = repmat([ones(Ny-1,1);0],Nx-1,1);
one_to_seven = repmat([one_to_seven;zeros(Ny,1)],Nz-1,1);

one_to_seven = one_to_seven(1:end-Ny-1);

two_to_eight = repmat([0;ones(Ny-1,1)],Nx-1,1);
two_to_eight = repmat([two_to_eight;zeros(Ny,1)],Nz-1,1);

two_to_eight = two_to_eight(1:end-Ny+1);

three_to_five = repmat([0;ones(Ny-1,1)],Nx-1,1);
three_to_five = repmat([zeros(Ny,1);three_to_five],Nz-1,1);
three_to_five = [three_to_five;zeros(Ny+1,1)];

four_to_six = repmat([ones(Ny-1,1);0],Nx-1,1);
four_to_six = repmat([zeros(Ny,1);four_to_six],Nz-1,1);
four_to_six = [four_to_six;zeros(Ny-1,1)];
%
% We can now combine all of this connections to form the adjacency matrix.
% We also need to keep track of the diagonals as well as which type of diagonal it is.
% There are two types of diagonals, small and large.
% The small diagonals are those that do not cross the center of the cube.
% The large diagonals do cross the center of the cube.
%

A = zeros(NP,NP);
small_diagonals = zeros(NP,NP);
large_diagonals = zeros(NP,NP);
if Nx>1
    A = A+diag(x_direction,Ny);
    if Ny>1
       A = A + diag(bl_to_tr_xy,Ny+1)+diag(tl_to_br_xy,Ny-1);
       small_diagonals = small_diagonals+...
                            diag(bl_to_tr_xy,Ny+1)+...
                            diag(tl_to_br_xy,Ny-1);
    end
    if Nz>1
       A = A+diag(tl_to_br_xz,Ny*Nx+Ny)+diag(bl_to_tr_xz,Nx*Ny-Ny);
       small_diagonals = small_diagonals+...
                        diag(tl_to_br_xz,Ny*Nx+Ny)+...
                        diag(bl_to_tr_xz,Nx*Ny-Ny);
    end
end
if Ny>1
    A = A+diag(y_direction,1);
    if z_direction>1
        A = A+diag(tr_to_bl_yz,Ny*Nx+1)+diag(tl_to_br_yz,Nx*Ny-1);
        small_diagonals = small_diagonals+...
                          diag(tr_to_bl_yz,Ny*Nx+1)+...
                          diag(tl_to_br_yz,Nx*Ny-1);
    end
end
if Nz>1
    A = A+diag(z_direction,Nx*Ny);
    if Nx>1 && Ny>1
       A = A+...
            diag(one_to_seven,Nx*Ny+Ny+1)+...
            diag(two_to_eight,Nx*Ny+Ny-1)+...
            diag(three_to_five,Nx*Ny-Ny-1)+...
            diag(four_to_six,Nx*Ny-Ny+1);
       large_diagonals = large_diagonals+...
            diag(one_to_seven,Nx*Ny+Ny+1)+...
            diag(two_to_eight,Nx*Ny+Ny-1)+...
            diag(three_to_five,Nx*Ny-Ny-1)+...
            diag(four_to_six,Nx*Ny-Ny+1);
    end
end
% The adjacency matricies should be symmetric we therefore 
% add the transposes of the upper triangular matrices to each other.
A = A+A';
small_diagonals = small_diagonals+small_diagonals';
large_diagonals = large_diagonals+large_diagonals';
end

