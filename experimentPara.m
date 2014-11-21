function expe = experimentPara()

%% don't change that

expe=struct();

%% change this

expe.numberOfFrames = 144;
expe.dt = 0.5;

expe.numberOfMovies = 2;
expe.indexOfFirstMovie = 21;

expe.numberOfColors = 2;
expe.colorNames = {'YFP','Cy32'};

imgDir = '/Volumes/Naef-Lab/Rosie/20141114_3t3_wt_or_3t3_Venus_transfection_plus_mKate_plus_or_3T3_Fucci/';

%% don't change that

mainDir =  mfilename('fullpath'); mainDir = fileparts(mainDir);
mainDir = [mainDir '/'];

expe.t = linspace(0,(expe.numberOfFrames-1) * expe.dt,expe.numberOfFrames);
N = expe.numberOfFrames;

