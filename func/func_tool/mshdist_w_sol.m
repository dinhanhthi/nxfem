function mshdist_w_sol(msh,phi,pathfile,name)
% Export phi to .sol file in order to use toolbox mshdist
% https://github.com/ISCDtoolbox/Mshdist
% CHECK: good
% Input: 
%   - mesh msh, 
%   - value of phi at vertices 
%   - path to export (example: /path/to/file/)
%   - name of file, without extension (example: phi)
% Output: file .sol

nP = size(msh.p,2); % number of vertices

pathfile = strcat(pathfile,name,'.sol');

fileID = fopen(pathfile,'w');
fprintf(fileID,'MeshVersionFormatted 2\n\n');
fprintf(fileID,'Dimension 2\n\n');
fprintf(fileID,'SolAtVertices\n%d\n',nP);
fprintf(fileID,'1 1\n\n');
fprintf(fileID,'%.15f\n',phi);
fprintf(fileID,'\nEnd');
fclose(fileID);

end