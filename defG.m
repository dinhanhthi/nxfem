function val=defG(phii,phij,typeG)
% Define function of phi_i, phi_j
% Input: phi_i, phi_j
% Output: 
    
if typeG==1
   val = phii.*phij; % g(X,Y)=XY
end

end