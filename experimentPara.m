function expe = experimentPara()

%% don't change that

expe=struct();

%% change this

expe.numberOfFrames = 200;
expe.dt = 0.5;

expe.numberOfMovies = 1;
expe.indexOfFirstMovie = 1;

expe.numberOfColors = 1;
expe.colorNames = {'FITC'};

expe.hasTrans = 0;
expe.transName = 'Trans';

expe.imgDir = '/Users/bieler/Desktop/moviesAndrea/onur/';

expe.mainDir = '/Users/bieler/Desktop/matlab/segmentation2';

%% don't change that

expe.t = linspace(0,(expe.numberOfFrames-1) * expe.dt,expe.numberOfFrames);
N = expe.numberOfFrames;



%%
