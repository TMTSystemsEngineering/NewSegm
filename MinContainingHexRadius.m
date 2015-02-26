function R = MinContainingHexRadius(points, tolerance)

% find max radius of points
maxR = 0;
for i=1:size(points,2)
    R = norm(points(:,i));
    maxR = max(maxR,R);
end;

% start at that radius and grow by tolerance until all points are inside
R = maxR - 3*tolerance;
R = R - mod(R,tolerance);
while ~AllPointsAreInside(R,points)
    R = R + tolerance;
end;
