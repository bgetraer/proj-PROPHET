experiment='RCP85';
expdir = ['/nobackupp18/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/' experiment];
rundir = fullfile(expdir,'runocean');
zipfilename=[experiment 'Coceanspinup_11_02_2024.zip'];
S = dir(rundir);

% dont include the boundary forcing files
S = S(~contains({S.name},'obs') & ~contains({S.name},'obw') & ~[S.isdir]);
filenames={S.name};
zip(fullfile(expdir,zipfilename),filenames,rundir);
