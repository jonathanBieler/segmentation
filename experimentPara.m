function expe = experimentPara()

%% don't change that

expe=struct();

%% change this

expe.numberOfFrames = 600;
expe.dt = 1/12;

expe.numberOfMovies = 1;
expe.indexOfFirstMovie = 1;

expe.numberOfColors = 2;
expe.colorNames = {'YFP','TexasRed'};

expe.hasTrans = 0;
expe.transName = 'Trans';

expe.imgDir = '/Volumes/Backup/2015_19_02_Rosie_Venus_H2B_Z_stacks_1/b4/';

expe.mainDir = '/Users/bieler/Desktop/matlab/histoneStuff/2015_19_02_Rosie_Venus_H2B_Z_stacks_1';

%% don't change that

expe.t = linspace(0,(expe.numberOfFrames-1) * expe.dt,expe.numberOfFrames);
N = expe.numberOfFrames;



%%
