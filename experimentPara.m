function expe = experimentPara()

%% don't change that

expe=struct();

%% change this

expe.numberOfFrames = 50;
expe.dt = 1/2;

expe.numberOfMovies = 1;
expe.indexOfFirstMovie = 51;

expe.numberOfColors = 2;
expe.colorNames = {'Cy33','FITC'};%{'YFP','TexasRed'};

expe.hasTrans = 1;
expe.transName = 'Trans'; %use Trans for bsf, see getOriginalImageName.m

expe.imgDir  = '/Users/bieler/Desktop/matlab/testMovies/cy3Fitch/';
expe.mainDir = '/Users/bieler/Desktop/matlab/segmentation3/';

%% don't change that

expe.t = linspace(0,(expe.numberOfFrames-1) * expe.dt,expe.numberOfFrames);
N = expe.numberOfFrames;



%%
