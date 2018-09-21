function mshdist_w_mesh(msh,pathfile,name)
% Export mesh data to .mesh file in order to use toolbox mshdist
% https://github.com/ISCDtoolbox/Mshdist
% CHECK: good
% Input: 
%   - mesh msh 
%   - path to export (example: /path/to/file/)
%   - name of file, without extension (example: phi)
% Output: file .mesh

nP = size(msh.p,2); % number of vertices
nT = size(msh.t,2); % number of triangles
points = [msh.p;zeros(1,nP)];
triangles = msh.t(1:3,:); % take 3 vertices
triangles = [triangles;zeros(1,nT)];

pathfile = strcat(pathfile,name,'.mesh');

fileID = fopen(pathfile,'w');
fprintf(fileID,'MeshVersionFormatted\n1\n\n');
fprintf(fileID,'Dimension\n2\n\n');
fprintf(fileID,'Vertices\n%d\n\n',nP);
fprintf(fileID,'%.15f %.15f %1d\n',points);
fprintf(fileID,'\nTriangles\n%d\n\n',nT);
fprintf(fileID,'%d %d %d %d\n',triangles);
fprintf(fileID,'\nEnd');
fclose(fileID);

end