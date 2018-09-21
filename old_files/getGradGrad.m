function gg = getGradGrad(tris,msh)
% This function is to find Grad*Grad of basis functions on triangles "tris"
% Input: - triangles "tris" (NCTs or CTs,...)
% Output: - vector contains values of GradGrad (3i x 3j x nTris)

nTris = size(tris,2); % number of triangles "tris"
% Gradient of function phi_i of the triangles (i is the vertex of triangle)
gradPhiTris = getGradPhi(tris,msh); % 2 coordinates x 3 vertices x nTris
gg = zeros(3,3,nTris);

for t=1:nTris
   for i=1:3
      for j=1:3
         gg(i,j,t) = dot(gradPhiTris(:,j,t),gradPhiTris(:,i,t));
      end
   end
end
    
end