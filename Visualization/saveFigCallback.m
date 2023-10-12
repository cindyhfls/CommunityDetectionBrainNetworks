function saveFigCallback()
    % Define your directory and filename
    saveDir = './tmp';  % Change this to your path
    fileName = 'myFigure.png';         % Change this to your filename
    fullPath = fullfile(saveDir, fileName);
    
    % Save the figure
    saveas(gcf, fullPath);
    disp(['Figure saved as ' fullPath]);
end
