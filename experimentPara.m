function expe = experimentPara()

%% don't change that

expe=struct();

%% change this

expe.numberOfFrames = 50;
expe.dt = 1/2;

expe.numberOfMovies = 1;
expe.indexOfFirstMovie = 51;

expe.numberOfColors = 1;
expe.colorNames = {'LUC'};%{'YFP','TexasRed'};

expe.hasTrans = 0;
expe.transName = 'Trans'; %use Trans for bsf, see getOriginalImageName.m

expe.imgDir  = '/Users/bieler/Desktop/matlab/testMovies/cy3Fitch/';
expe.mainDir = '/Users/bieler/Desktop/matlab/segmentation3/';

%% don't change that

expe.t = linspace(0,(expe.numberOfFrames-1) * expe.dt,expe.numberOfFrames);
N = expe.numberOfFrames;



%%
