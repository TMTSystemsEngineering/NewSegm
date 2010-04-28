function radius = circumradius(center,points);
% function radius = circumradius(center,points);
% E. Ponslet, 2007, elp@ericandlucie.com
% Calculates radius of enclosing circle of set of points, given its center
% 
% INPUTS:
%    * center(2 x 1): in plane coordinates of a candidate center for the enclosing circle
%    * points(2 x npoints): in plane coordinates of any number of points to be enclosed
% OUTPUTS:
%   * radius (scalar): radius of the circle centered at 'center' that encloses all 'points'

dxy=points-center*ones(1,size(points,2));  % calculate vectors from 'center' to all 'points'
radius=max(sqrt(dxy(1,:).^2+dxy(2,:).^2));  % calculate radius of circumscribed circle as length of the longest of those vectors
