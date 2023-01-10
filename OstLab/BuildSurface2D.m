function circle_surface = BuildSurface2D(Nx,Ny,R,dist,n_dims)
%BUILDFLOOR Build 2D floor consisting of randomly distributed circles.
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
% if length(R)<N
%     R = repmat(R,N,1);
%     R = R(1:N,1); 
% end
circle_centers = zeros(Nx,Ny,n_dims);
[x,y]=meshgrid(0:R:(Nx-1)*R,0:R:(Ny-1)*R);
circle_centers(:,:,1)=x';
circle_centers(:,:,2)=y';
circle_centers(2:end,2:end,1:2) = circle_centers(2:end,2:end,1:2)+cumsum(dist*randn(Nx-1,Ny-1,2).*R,1);
% keyboard
% Center the first circle on the origin.
circle_centers = reshape(circle_centers, [Nx*Ny,n_dims]);
circle_surface = [circle_centers repmat(R,Nx*Ny,1)];

% Last dimension of the circle_centers



end