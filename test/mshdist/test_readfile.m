% fileID = fopen('test2.sol','r');
% nV = 6;
% FormatString = repmat('%d',1,nV);  
% C = textscan(fileID,FormatString,'Delimiter','\n','MultipleDelimsAsOne',1,'HeaderLines',6);
% fclose(fileID);
% disp(C);

b=5;
fileID = fopen('test.sol','r');
% nV = textscan(fileID,'%d',1,'Delimiter','\n','HeaderLines',7);
% nV = nV{1};
C = textscan(fileID,'%f',10,'Delimiter','\n','HeaderLines',10);
fclose(fileID);
A=C{1}
