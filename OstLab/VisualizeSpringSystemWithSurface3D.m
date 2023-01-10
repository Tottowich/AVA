function  VisualizeSpringSystemWithSurface3D(X,A,circle_surface,record)
%VISUALIZECHEESE Visualize the simulated spring system'cheese'
%   Function used to animated the system of springs.
%
% Author: Theodor Jonsson
% Date: 05/01/2023
% INPUT
%
%   X - (mat) matrix of size (t_steps, NP, n_dims)
%
%   A - (mat) Adjacency matrix of the spring system.
%
%   circle_surface - (mat) Matrix of shape (N x n_dims+1). For each circle
%                          there is center (x,y,z) + radius.
%
[t_steps,NP,n_dims] = size(X);
set(gcf,'Position',get(0,'Screensize')); % Maximize the window for quality
figure(1)
if record
    % Initialize animation file before main loop
    MOVE = VideoWriter('ExampleVideo.avi'); % Create an video struct, "MOVE.",
    % with the output file ExampleVideo.avi
    MOVE.Quality = 50; % Set the quality (0-100)
    MOVE.FrameRate = 25; % Set frames per seond to PAL standard (animation speed)
    open(MOVE);
end
set(gca,'nextplot','replacechildren');

% Set the limits to avoid window stuttering.
minminx = min(X(:,:,1),[],'all');
minminy = min(X(:,:,2),[],'all');
minminz = min(X(:,:,3),[],'all');
maxmaxx = max(X(:,:,1),[],'all');
maxmaxy = max(X(:,:,2),[],'all');
maxmaxz = max(X(:,:,3),[],'all');

xlim([minminx-1,maxmaxx+1])
ylim([minminy/2-1,maxmaxy+1])
zlim([minminz/2-1,maxmaxz+1])

daspect([1,1,1]);
grid on
hold on
% Initialize the grid of springs.
x = squeeze(X(1,:,1));
y = squeeze(X(1,:,2));
z = squeeze(X(1,:,3));
springs_linespec = "k.-";
SPRINGS = plot(graph(A),springs_linespec,'XData',x,'YData',y,'ZData',z,'NodeLabel',{});
view(90,0)
% Display the circular floor.
N_cirlces = length(circle_surface);
[sx,sy,sz] = sphere(40); % Sphere for plotting
for n = 1:N_cirlces
    % Create a cirlce as a rectangle.
    % Since the rectangle has position based bottom left, and top right.
    % circle_surface contains the center of the circle, find ["bottom
    % left","top right"] for the circle.
    pos = circle_surface(n,1:n_dims);
    r = circle_surface(n,end);
    spheres(n) = surf(sx*r+pos(1),sy*r+pos(2),sz*r+pos(3),'EdgeColor','k','EdgeAlpha',0.05,'FaceColor',[0,0.8,0.35],'FaceAlpha',0.4);
%     pos = [c(1)-r,c(2)-r,2*r,2*r];
%     spheres(n) = rectangle(Curvature=[1,1],Position=pos);
    % Create the rectangle which will be a circle.

end
% Begin main time loop
for t = 2:size(X,1)
    x = squeeze(X(t,:,1));
    y = squeeze(X(t,:,2));
    z = squeeze(X(t,:,3));
    % Update the grid.
    set(SPRINGS,'XData',x,'YData',y,'ZData',z)
    if record
        writeVideo(MOVE,getframe(gcf)); % Get a snapshot of the active figure frame
    else
        pause(20/length(X)); % The video takes roughly 20 seconds
    end
    % End of main time loop
end
hold off
if record
    close(MOVE); % Close and save the avi-file
end

end % Function

