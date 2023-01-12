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
                                               % (for horizontal connections)

horizontal = horizontal(1:end-1);              % Remove the last value

tl_to_br = [0; horizontal(1:(Nc*(Nr-1)))]; % Make the second diagonal
                                   % vector (for top left to bottom right)

vertical = ones(Nc*(Nr-1), 1);         % Make the third diagonal vector
                                       %   (for vertical connections)

bl_to_tr = tl_to_br(2:end-1); % Make the fourth diagonal vector
                                       %   (for bottom left to top right)
A = diag(horizontal, 1)+...            % Add the diagonals to a zero matrix
    diag(vertical, Nc)+...
    diag(bl_to_tr, Nc+1)+...
    diag(tl_to_br, Nc-1);
%A = A-diag(horizontal,1);
A = A+A.'; % This makes the Adjacency matrix symmetric.
% Also return the diagonals to make the diagonal springs longer.
diagonals = diag(tl_to_br, Nc-1)+diag(bl_to_tr, Nc+1);
diagonals = diagonals+diagonals';
end

