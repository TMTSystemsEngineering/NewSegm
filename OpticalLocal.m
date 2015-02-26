function z=OpticalLocal(x,y,k,K,Center,Onex,Oney,Onez);
% function z=OpticalLocal(x,y,k,K,Center,Onex,Oney,Onez);
% Eric Ponslet, 2007, elp@ericandlucie.com
% Returns the local z coordinate of a point in the M1 optical surface, given its local x,y coordinates
%  * the local frame is defined by Center (3x1), and Onex, Oney, Onez (3x1), the unit vectors of the local frame 
%  * the optical surface is defined by K and K, according to the equation (in global coordinates) Z= {k-sqrt[k^2-(K+1)*R^2]} / (K+1)
% Inputs x and y can be arrays (to compute several points at once), but they must be the same size; z will be same size as x and y.
% Technique: solve second order equation for local z coordinate that results from subsituting the equation 
%            of a line parallel to local z originating from x,y, into the equation of the optical surface

un=ones(size(x));

% calculate intermediate terms
DX = Center(1)*un + Onex(1)*x + Oney(1)*y;
DY = Center(2)*un + Onex(2)*x + Oney(2)*y;
DZ = Center(3)*un + Onex(3)*x + Oney(3)*y;

% calculate coefficients of 2nd order polynomial
a= un*( (K+1)*Onez(3)^2 + Onez(1)^2 + Onez(2)^2 );
b= 2*DX*Onez(1) + 2*DY*Onez(2) + 2*DZ*Onez(3)*(K+1) - 2*k*Onez(3)*un;
c= DX.^2 + DY.^2 + (K+1)*DZ.^2 - 2*k*DZ;

% calculate root
z=(-b-sqrt(b.^2-4.*a.*c))./(2*a);
