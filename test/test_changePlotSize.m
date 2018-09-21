%% cf. http://bit.ly/2OK6TUT


% %% depend on the size on screen & screen resolution
% surf(peaks)
% set(gcf,'PaperPositionMode','auto')
% print('test/savePlot/PeaksSurface','-dpng','-r0')


%% with specific dimension size
% get the default value
oldpaperunits = get(gcf,'PaperUnits');
oldpaperpos = get(gcf,'PaperPosition');

% preserve axis and limits
ax = gca; 
ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';

surf(peaks);
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 8 6];
print('test/savePlot/PeaksSurface_PaperUnits','-dpng','-r0');

fig.PaperUnits = oldpaperunits;
fig.PaperPosition = oldpaperpos;