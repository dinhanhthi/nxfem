fileName='thithi.txt';
fileID = fopen(fileName,'w');
        fprintf(fileID,'%s,\n','aaa');
% fclose(fileID); 
a=5; b=6;
c=a+b;
fprintf(fileID,'%s,\n','bbb');
fclose(fileID); 