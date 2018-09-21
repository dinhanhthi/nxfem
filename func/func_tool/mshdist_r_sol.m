function phinew = mshdist_r_sol(phi,pathfile,name)
% read the file .sol
% CHECK: good
% Input:
%   - the old phi, 
%   - path to export (example: /path/to/file/)
%   - name of file, without extension (example: phi)
% Output: new phi (1 x number of vertices)

nV = size(phi,2);

pathfile = strcat(pathfile,name,'.sol');

fileID = fopen(pathfile,'r');
C = textscan(fileID,'%f',nV,'Delimiter','\n','HeaderLines',8);
fclose(fileID);

phinew = C{1};
phinew = phinew'; % convert to row array

end