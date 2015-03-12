function A = print2array(fig, res)
% Generate default input arguments, if needed
if nargin < 2
    res = 300;
    if nargin < 1
        fig = gcf;
    end
end
% Set paper size
set(fig, 'PaperPositionMode', 'auto');
% Generate temporary file name
tmp_nam = ['temp_image.tif'];
% Print to tiff file


print('-noui' , ['-f' num2str(fig)], ['-r' num2str(res)], '-dtiff', tmp_nam);

% Read in the printed file
A = imread(tmp_nam);
% Delete the file
delete(tmp_nam);
return