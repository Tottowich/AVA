
% Function to plot 3D-function
% Should take as input:
% f_X - anonymous function of y and z
% f_Y - anonymous function of x and z
% f_Z - anonymous function of x and y
% x - vector of x-values
% y - vector of y-values
% z - vector of z-values
% Should return:
% h - handle to the plot
function [h,X,Y,Z] = plot3D(f_X,f_Y,f_Z,x,y,z)
    % Create a meshgrid
    [x,y] = meshgrid(x);
    z = y;

    %size(z)
    % Calculate the values of the functions
    X = f_X(x,y,z);
    Y = f_Y(x,y,z);
    Z = f_Z(x,y,z);
    % Plot the function
    h = surf(X,Y,Z,"DisplayName","f(x,y,z)");
    % Label the axes
    xlabel('x');
    ylabel('y');
    zlabel('z');
    % Add lines through the origin
    hold on;
    plot3([0 0],[0 0],[min(min(Z)) max(max(Z))],'g',"DisplayName","z-axis","LineWidth",2);
    plot3([0 0],[min(min(Y)) max(max(Y))],[0 0],'r',"DisplayName","y-axis","LineWidth",2);
    plot3([min(min(X)) max(max(X))],[0 0],[0 0],'b',"DisplayName","x-axis","LineWidth",2);
    scatter3(0,0,0,'r',"DisplayName","Origin","LineWidth",4);
    legend()
    % Create normals of the surface:
    % The normal is the vector perpendicular to the surface
    % The normal is the cross product of the partial derivatives

    % Calculate the partial derivatives
    dX_dx = diff(X,1,1);
    dX_dy = diff(X,1,2);
    dY_dx = diff(Y,1,1);
    dY_dy = diff(Y,1,2);
    dZ_dx = diff(Z,1,1);
    dZ_dy = diff(Z,1,2);

    dX_dx = [dX_dx; dX_dx(end,:)];
    dX_dy = [dX_dy dX_dy(:,end)];
    dY_dx = [dY_dx; dY_dx(end,:)];
    dY_dy = [dY_dy dY_dy(:,end)];
    dZ_dx = [dZ_dx; dZ_dx(end,:)];
    dZ_dy = [dZ_dy dZ_dy(:,end)];
    % Calculate the normals
    N_x = dY_dy.*dZ_dx - dY_dx.*dZ_dy;
    N_y = dX_dx.*dZ_dy - dX_dy.*dZ_dx;
    N_z = dX_dy.*dY_dx - dX_dx.*dY_dy;
    % Normalize the normals
    % Calculate the areas of the rectangles
    A = sqrt(N_x.^2 + N_y.^2 + N_z.^2);
    % Calculate the area of the surface
    area = sum(sum(A));
    fprintf("Area of the surface: %f",area)
    % Normalize the normals
    N = sqrt(N_x.^2 + N_y.^2 + N_z.^2);
    N_x = N_x./N;
    N_y = N_y./N;
    N_z = N_z./N;
    % Plot the normals with length 10
    % Increase length of quiver arrows
    % Only plot every kth arrow
    k = 15;
    quiver3(X(1:k:end,1:k:end),Y(1:k:end,1:k:end),Z(1:k:end,1:k:end),N_x(1:k:end,1:k:end),N_y(1:k:end,1:k:end),N_z(1:k:end,1:k:end),1,"DisplayName","Normal vector");
    hold off;
    % Set the aspect ratio
    %
    daspect([1 1 1]);
    % Fix axis limits
    axis manual;
    % Set axis to equal range and add some margin
    % max_range = max([max(max(X)) max(max(Y)) max(max(Z))])+1;
    % min_range = min([min(min(X)) min(min(Y)) min(min(Z))])-1;
    % axis([min_range max_range min_range max_range min_range max_range]);
    % Zoom out a bit
    camzoom(1.5);
    % Bring the plot to the front
    uistack(h,'top');
    figure(2);
    % Plot X, Y, and Z as images (2D plots)
    subplot(1,3,1);
    imagesc(X);
    title('X');
    subplot(1,3,2);
    imagesc(Y);
    title('Y');
    subplot(1,3,3);
    imagesc(Z);
    title('Z');

end