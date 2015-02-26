function [RH_radius,RH_clocking,Irregularity]=FitHexNew(Vertices,a)
% function [RH_radius,RH_clocking,Irregularity]=FitHex(Vertices,a)
% E. Ponslet, February 2010, elp@ericandlucie.com
% Fits a regular hexagon to a set of six vertices 'Vertices(2,6)' in the least square sense.
% The best fit regular hexagon (BFRH) is allowed to scale and clock, but not recenter.
% Returns the radius and clocking angle of the best fit hexagon, and the RSS vertex position error, a measure of irregularity
%
% INPUTS:
%   * Vertices (array, 2*6): X and Y coordinates of the 6 vertices of the irregular hexagon
%   * a: starting guess for the radius of the BFRH
% OUTPUTS:
%   * RH_radius (scalar): radius of the BFRH
%   * RH_clocking (scalar): clocking angle of the BFRH (angle from X axis to vertex #1 of the BFRH), in radians
%   * Irregularity (scalar): RSS of the distances from each vertex of the BFRH to the corresponding vertex of the irregular hexagon defined by 'vertices' 

options = optimset('Display','off','TolFun',sqrt(6)*1e-14);  % set options for optimization routine, and termination tolerance

% minimize the square root of the sum of the squares of distances between actual vertices and best fit vertices, measured in the local XY plane
x0=[a;0];      % starting guess
SEGvertices=[Vertices(1,:)';Vertices(2,:)'];  % rearrange 2D vertex coordinates into long (12x1) column vector (to ease calculations in HexNew)
[x,resnorm]=fminsearch(@HexNew,x0,options,SEGvertices);  % NL fit routine, returns x=[radius,clocking] and Sum Square residual; passes "Vertices" to evaluation function HexNew
RH_radius=x(1);
RH_clocking=x(2);

% calculate RMS irregularity
Irregularity=resnorm/sqrt(6);  % this is an RMS vertex position error

return
