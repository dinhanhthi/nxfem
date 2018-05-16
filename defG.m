function val=defG(u,typeG)
% Define function g(u) and its derivative
% Input: u vector of scalar
% Output: value of defG depending on typeG

% Must change both of gu and dgu!!!!
gu = u.^2; % change the form here
dgu = 2*u; % g'(u)

switch typeG
    case 1 % g(u)
        val = gu;
    case 2 % g=0
        val = 0.; % just for testing
    case 3 % g'(u)
        % newton method for finding u (linda's test case,cf.getGMuNewton.m)
        val = dgu;
    case 4 % u*g'(u)-g(u)
        % newton method for finding u (linda's test case,cf.getGMuNewton.m)
        % cf. int_Omg g(uold)*phi*phi in .getGMuNewton.m
        val = u.*dgu-gu;
    case 5 % u
        val = 1.; % for Chopp07
end

end