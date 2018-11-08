x=[1,2,3];
y=[1,5,1];
y2=[-3,-1,0];

% oldpaperunits = get(gcf,'PaperUnits');
% oldpaperpos = get(gcf,'PaperPosition');
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 8 6];

% f=figure; 
% set(f, 'Visible', 'off');
% plot(x,y,'-r');
% fileName = 'test.png';
% % oldpaperunits = get(gcf,'PaperUnits');
% % oldpaperpos = get(gcf,'PaperPosition');
% % fig = gcf;
% % fig.PaperUnits = 'inches';
% % fig.PaperPosition = [0 0 8 6];
% print(fileName,'-dpng','-r0');
% % hold on
% 
% plot(x,y2,'-b');
% fileName = 'test.png';
% print(fileName,'-dpng','-r0');
% 
% g=figure; 
% set(g, 'Visible', 'off');
% plot(x,y,'-r');
% fileNameg = 'testg.png';
% print(fileNameg,'-dpng','-r0');
% close(g);

%% try with saveas

f = figure('visible', 'off');

if ~isempty(f)
    disp('f exists');
else
    disp('f doesnt exist');
end

% plot(x,y,'-r');
% oldpaperunits = get(gcf,'PaperUnits');
% oldpaperpos = get(gcf,'PaperPosition');
% fig = gcf;
% fig.PaperUnits = 'inches';
% fig.PaperPosition = [0 0 8 6];
% saveas(f,'fr','png');

close(f);

