%% Testing Section
main()


%% Main
function main()

end


%% Molecule Collision Detection and Correction:

function collisions = detect_element_collision(data, size)
% Returns an upper triangular logical matrix where the row and column
% denote the identity of the two particles, and the value is true if these
% two particles have collides.
% Takes:
%   data: A matrix of 2D position vectors where each column is a 
%       different point
%   size: The size of each molecule
% Returns:
%   collisions: A logical upper triangular matrix of identified collisions

    collisions = boolean(zeros(width(data))); % Logical array
    size_sqr = size^2;
    
    for j = 1:(width(data) - 1) % Do not check the last piece of data
    
        dr = data(:, (j + 1):end ) - data(:,j);
        L_sqr = sum(dr .* dr);
        collisions(j,(j+1):end) = L_sqr <= size_sqr;
    
    end
    
end
function [xP1, xP2] = correct_element_position(x1, v1, x2, v2, size)
% If the particles significantly overlap, then this function will move the
% particles backwards in time until they no longer overlap.
% Takes:
%   x1, x2: The position vectors
%   v1, v2: The velocity vectors
%   size: The radius of the particle
% Returns:
%   xP1, xP2: The corrected position vectors of the particles

    % Calculate Necessary Variables:
    dr = x1 - x2;
    dv = v1 - v2;
    R_sqr = size^2;
    
    % dt: Can be solved for using a quadratic equation
    dt = roots([dot(dv,dv),-2 * dot(dr,dv), dot(dr,dr) - (4 * R_sqr)]);
    
    % Choose the positive root:
    if dt(1) > 0
        dt = dt(1);
    else
        dt = dt(2);
    end
    
    % Make Correction:
    xP1 = x1 - (dt * v1);
    xP2 = x2 - (dt * v2);

end