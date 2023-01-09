function [A,diagonals] = GridAdjacencyMatrix3D(Nx,Ny,Nz)
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
horizontal = repmat([ones(Nx-1, 1); 0], Ny, 1);% Make the first diagonal vector
                                               % (for horizontal connections)

horizontal = horizontal(1:end-1);              % Remove the last value

tl_to_br = [0; horizontal(1:(Nx*(Ny-1)))]; % Make the second diagonal
                                   % vector (for top left to bottom right)

vertical = ones(Nx*(Ny-1), 1);         % Make the third diagonal vector
                                       %   (for vertical connections)

bl_to_tr = tl_to_br(2:end-1); % Make the fourth diagonal vector
                                       %   (for bottom left to top right)
A = diag(horizontal, 1)+...            % Add the diagonals to a zero matrix
    2*diag(tl_to_br, Nx-1)+...
    3*diag(vertical, Nx)+...
    4*diag(bl_to_tr, Nx+1);
A = A+A.'; % This makes the Adjacency matrix symmetric.
diagonals = diag(tl_to_br, Nx-1)+diag(bl_to_tr, Nx+1);
diagonals = diagonals+diagonals';
% figure(4)
% c = 256/4;
% image(A*c)
% % 3D Should contain 13 different diagonals of the matrix.
% % This is because it can be a total of 26 connections to a single particle.
% x_direction = null(1);
% y_direction = repmat([ones(Nx-1, 1); 0], Ny*Nz, 1);
% % Remove last
% y_direction = y_direction(1:end-1);
% z_direction = null(1);
% keyboard
end

