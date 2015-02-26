function F=HexNew(x,Vertices)
% function F=HexNew(x,Vertices)
% E. Ponslet, February 2010, elp@ericandlucie.com
% Calculates RSS of 2D distances between Vertices of an irregular hexagon and those of a regular hexagon of radius x(1) and clocking angle x(2)
%
% INPUTS
%   * x (2*1): contains the radius (x(1)) and clocking angle (x(2), in radians) of the regular hexagon
%   * Vertices (12*1): coordinate of the vertices of the irregular hexagon, rearranged into a colum vector [X1;X2;...X6;Y1;Y2;...Y6]
% OUTPUTS
%   * F (scalar): RSS of distances from vertices of irregular hexagon to corresponding vertices of regular hexagon

Radius=x(1);
Rotation=x(2);

% vertex angles about center for a regular hexagon rotated by angle x(2) 
angles=[0:5]'*pi/3+Rotation;
% Coordinates of vertices of regular hexagon, arranged into one 12x1 column vector
BFvertices=[Radius*cos(angles);Radius*sin(angles)];
% calculate RSS of distances from Vertices to BFvertices
F=norm(BFvertices-Vertices);
