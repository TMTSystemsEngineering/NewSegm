function [Z,R,dZdR]=OpticalSurface(X,Y,k,K)
% function [Z,R,dZdR]=OpticalSurface(X,Y,k,K)
% Eric Ponslet, 2007
% Calculates the global Z coordinate of point(s) in the optical surface at global coordinates X, Y.
% Also returns the Radial coordinate of point(s) X,Y, and the slope dZ/dR of the optical surface at point X,Y,Z 
%
% INPUTS:
%   * X,Y (scalars or arrays): coordinates of a point (or several points) in M1CRS; If X and Y are arrays, they must have the same dimensions.
%   * k (scalar): paraxial radius of curvature of the optical surface
%   * K (scalar): conic constant of the optical surface
% OUTPUTS:
%   * Z (scalar or array, same size as X & Y): Z coordinate in M1CRS of point(s) in the optical surface defined by X and Y
%   * R (scalar or array, same size as X & Y): R coordinate in M1CRS of point(s) defined by X and Y
%   * dZdR (scalar or array, same size as X & Y): radial slope of the optical surface at point(s) X,Y,Z

R=sqrt(X.^2+Y.^2);
unit=ones(size(R));

% altitude at that point
Z=( k*unit - sqrt(k^2*unit-(K+1)*R.^2) ) / (K+1);  % full equation of M1 surface

% local slope in radial direction
dZdR = R ./ sqrt(k^2*unit-(K+1)*R.^2);   % derivative dZ/dR
