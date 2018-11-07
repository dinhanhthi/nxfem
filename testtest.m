x=[1,2,3];
y=[1,5,1];
y2=[-3,-1,0];

f=figure; 
set(f, 'Visible', 'off');
plot(x,y,'-r');
fileName = 'test.png';
print(fileName,'-dpng','-r0');
hold on

plot(x,y2,'-b');
fileName = 'test.png';
print(fileName,'-dpng','-r0');

g=figure; 
set(g, 'Visible', 'off');
plot(x,y,'-r');
fileNameg = 'testg.png';
print(fileNameg,'-dpng','-r0');
close(g);

close(f);

