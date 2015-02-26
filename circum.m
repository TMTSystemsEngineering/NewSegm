function [center,radius]=circum(points);
% function [center,radius]=circum(points);
% Eric Ponslet, March 2010, elp@ericandlucie.com
% Finds the center and radius of a minimum radius circle circumscribed to a set of points in 2D space
% Technique: use numerical optimization routine to find the x and y coordinates of a the center of a circle that minimizes its radius while enclosing all the 'points'
% 
% INPUTS:
%   * points (2*n): x and y coordinates of the points to be enclosed 
% OUTPUTS:
%   * center (2*1): x and y coordinates of the center of the minimum circumscribed circle
%   * radius (scalar): radius of the minimum circumscribed circle

%Initial guess for center
center=[0;0];

% iterative refinement
options=[];  % all default optionms for minimization algorithm
[center,radius] = fminsearch(@circumradius,center,options,points);  % unconstrained minimization
