% Program to test the function plot3D.m
clear all
close all
R = 1;
a = 2;
b = 3;
c = 4;
height = 4;
% f_X = @(x,z,y) a*R.*sin(x).*cos(y/2);
% f_Y = @(x,z,y) b*R.*sin(x).*sin(y/2);
% f_Z = @(x,z,y) c*R.*cos(x);
f_X = @(x,z,y) x;
f_Y = @(x,z,y) cos(z);
f_Z = @(x,z,y) y-sin(x);
x = linspace(0,2*pi,1000);
y = linspace(0,2*pi,200);
z = linspace(0,2*pi,200);

[h,X,Y,Z] = plot3D(f_X,f_Y,f_Z,x,y,z);

% Add git lfs tracking to this file:
% git lfs track "plot3D.m"