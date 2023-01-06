function [A,diagonals] = GridAdjacencyMatrix(Nr,Nc)
%GRIDADJACENCYMATRIX adjacency matrix of all 8 neighbors in 2D-grid.
%
% NOTE THE CODE IS INSPIRED BY 
%
% https://stackoverflow.com/questions/3277541/construct-adjacency-matrix-in-matlab/3283732#3283732
% 
% Create the adjacency matrix with all neighbors connected.
%
% Example of such grid:
%
%   +---+---+
%   | X | X |
%   +---+---+
%   | X | X |
%   +---+---+
% Plusses are particles lines, dashes and X's are connections.
%
% INPUT
%
%   Nr - (int) Number of rows in the grid.
%
%   Nc - (int) Number of columns in the grid.
%
% OUTPUT
%
%   A - (mat) The adjacency matrix of all connected points.
%
%   diagonals - (mat) The points of diagonal springs.
%
horizontal = repmat([ones(Nc-1, 1); 0], Nr, 1);% Make the first diagonal vector
                                             %   (for horizontal connections)

horizontal = horizontal(1:end-1);                % Remove the last value

righ_to_left = [0; horizontal(1:(Nc*(Nr-1)))]; % Make the second diagonal vector
                                             %   (for anti-diagonal connections)

vertical = ones(Nc*(Nr-1), 1);                 % Make the third diagonal vector
                                             %   (for vertical connections)

left_to_right = righ_to_left(2:end-1);            % Make the fourth diagonal vector
                                             %   (for diagonal connections)
A = diag(horizontal, 1)+...                  % Add the diagonals to a zero matrix
      diag(righ_to_left, Nc-1)+...
      diag(vertical, Nc)+...
      diag(left_to_right, Nc+1);
A = A+A.'; % This makes the Adjacency matrix symmetric.
diagonals = diag(righ_to_left, Nc-1)+diag(left_to_right, Nc+1);
diagonals = diagonals+diagonals';
end

