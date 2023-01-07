function circle_surface = BuildSurface(N,R,dist,n_dims)
%BUILDFLOOR Build floor consisting of randomly distributed circles.
%
% Creates a matrix describing the floor made of circles.
% The circles will be distributed somewhat randomly (dist*R*rand) 
% between each center, rand ranges from 0-1. The structure of the matrix is
% as follows;
%   size(circle_surface) == [N,n_dims+1], for each circle there is a center
%   position + radius => 
%       circle_surface = (circle number, [x,y,(z if3D),r])
%   
% INPUT
%
%   N - (int) Number of circles to create the floor of.
%
%   R - (int/vec) Radius for each circle.
%
%   dist - (float) Maximum proportion of R to seperate the center of the
%                  circles with.
%   n_dims - (int) Number of dimensions in the system.
%   
% OUTPUT
%
%   circle_surface - (mat) Matrix of shape (N x n_dims+1). For each circle
%                          there is center (x,y,z) + radius.
%
if length(R)<N
    R = repmat(R,N,1);
    R = R(1:N,1); 
end
circle_centers = zeros(N,n_dims);
circle_centers(2:end,1) = cumsum(R(2:end)+dist*randn(N-1,1).*R(2:end),1);
% Center the first circle on the origin.
circle_surface = [circle_centers R];

% Last dimension of the circle_centers



end