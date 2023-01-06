function circle_surface = BuildSurface(N,R,dist,n_dims)
%BUILDFLOOR Build floor consisting of randomly distributed circles.
%
% Creates an array of circle structures. The circles will be distributed
% somewhat randomly (dist*R*rand) between each center, rand ranges from
% 0-1.
% Each circle has the following properties.
if length(R)<N
    R = repmat(R,N,1);
end
circle_centers = zeros(N,n_dims);
circle_centers(2:end,1) = cumsum(R(1:N-1)+dist*randn(N-1,1).*R(1:N-1),1);
% Center the first circle on the origin.
circle_centers = circle_centers-ones(1,n_dims)*R(1)/2;;
for n = 1:N
    circle_surface(n).center = circle_centers(n,:);
    circle_surface(n).radius = R(n);
end


end