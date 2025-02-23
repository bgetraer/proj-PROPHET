% Run minimum working example of pickup with no coupling
experimentname='Paris2C';
niter0=933120;
% Setup directory from input, build, and climate forcings {{{
disp('Initialize rundir:');
builddir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/MITgcm_initialization/build';
inputdir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/MITgcm_initialization/input';
runoceandir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/Paris2C/runocean';
mwedir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/MWE';
rundir='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/MWE/run';
if isdir(rundir),rmdir(rundir,'s'); end
if ~isdir(mwedir),mkdir(mwedir); end
mkdir(rundir);
cd(rundir);

mitfilepath='/nobackup/bgetraer/issmjpl/proj-getraer/proj-PROPHET/experiments/MITgcm_initialization/Models/PROPHET_mitgcm_init_RuntimeOptionsCoupled.mat';
load(mitfilepath); % mit structure

% RUNOCEANDIR FILES
disp('   copying pickup files');
pickupSuffrolling='ckptB';
pickupSuffpermanent=sprintf('%010i',niter0);
source=fullfile(runoceandir,['pickup.' pickupSuffrolling '.meta']);
destination=['pickup.' pickupSuffpermanent '.meta'];
copyfile(source,destination);
disp(['         - ',destination]);
source=fullfile(runoceandir,['pickup.' pickupSuffrolling '.data']);
destination=['pickup.' pickupSuffpermanent '.data'];
copyfile(source,destination);
disp(['         - ',destination]);

% INPUT FILES
filelist={'eedata','data','data.cal','data.diagnostics','data.exf','data.pkg','data.shelfice',...
	'bathy.bin','delr.bin','draft.bin','sref.bin','tref.bin'};
ln_filelist(inputdir,filelist,rundir); % link all files from inputdir to runcoupleddir
% renamed link to data.obcs
disp('   linking OBCS file');
link='data.obcs';
target=fullfile(inputdir,['data.obcs' experimentname]);
system(['ln -s ' target ' ' link ]);
disp(['         - ',link]);

% OBCS FILES
S=dir(mit.forcing.Ddir); % forcing files
S=S(contains({S.name},experimentname) & (contains({S.name},'2012') | contains({S.name},'2013'))); % the only ones we need
filelist={S.name};
ln_filelist(mit.forcing.Ddir,filelist,rundir); % link them

% link mitgcmuv executable
disp('   linking MITgcm executable');
link='mitgcmuv';
target=fullfile(builddir,'mitgcmuv');
system(['ln -s ' target ' ' link ]);
disp(['         - ',link]);
% }}}
% Set transient options {{{
% timestepping options
disp(['Set runtime options in data file']);
% Run Start and Duration
mit.inputdata.PARM{3}.startTime = mit.timestepping.spinupduration;                                % run start time for this integration (s)
mit.inputdata.PARM{3}.nIter0    = 0;                                % starting timestep iteration number
mit.inputdata.PARM{3}.nTimeSteps= 1; % number of timesteps to execute
% Restart/Pickup Files
mit.inputdata.PARM{3}.pChkptFreq=mit.inputdata.PARM{3}.deltaT; % leave pickup file at each timestep

datafilepath=fullfile(rundir,'data');
write_datafile(datafilepath, mit.inputdata.PARM, 'MODEL PARAMETERS');
% diagnostics options
disp(['Set runtime options in data.diagnostics file']);
mit.inputdata.DIAG{1}.N(1).frequency=0;
mit.inputdata.DIAG{1}.N(2).frequency=0;
mit.inputdata.DIAG{1}.N(3).frequency=1;
datafilepath=fullfile(rundir,'data.diagnostics');
write_datafile(datafilepath, mit.inputdata.DIAG, 'DIAGNOSTICS RUNTIME PARAMETERS');
% }}}


% run the MITgcm executable with MPI
command = ['mpirun -np 140 ./mitgcmuv > out 2> err'];
system(command);

function ln_filelist(parentdir,filelist,targetdir) % {{{
% LN_FILELIST soft-links a list of files located in parentdir to targetdir
   if ~isdir(parentdir)
      error('parentdir must be a directory');
   elseif any(~isfile(fullfile(parentdir,filelist)))
      error('filelist contains files which do not exist in parentdir');
   elseif  ~isdir(targetdir)
      error('targetdir must be a directory');
   end
   % link the files
   disp(['   linking ' num2str(numel(filelist)) ' files from ']);
   disp(['       ' parentdir ' to']);
   disp(['       ' targetdir]);
   for i=1:numel(filelist)
      file_path=fullfile(parentdir,filelist{i}); % file location
      command = ['ln -s ' file_path ' ' targetdir];
      if numel(filelist)<20
         disp(['         - ' filelist{i}]);
      elseif i==1
         disp(['         ...']);
      end
         system(command);
   end
end   % }}}
function write_datafile(fname,C,head) % {{{
   % WRITE_DATAFILE writes structures in C to fname {{{
   % INPUT: fname   string file to write
   %        C       cell array of structures (P) to write
   %        head    string of file header
   %
   % OUTPUT: writes to file with the following form:
   %
   % # head
   % # C{1}.description
   %  &C{1}.header
   %  C{1}.field1=value1,
   %  C{1}.field2=value2,
   %  ...
   %  &
   %
   % # C{2}.description
   %  &C{2}.header
   %  C{2}.field1=value1,
   %  C{2}.field2=value2,
   %  ...
   %  &
   %  ...
   % }}}
   disp(['    writing namelist file to ' fname])
   fileID = fopen(fname,'w');
   fprintf(fileID,'# %s\n',head); % write descriptive file header

	% loop through structures in C
   for i=1:length(C)
      % write a descriptive comment if it exists {{{
      if isfield(C{i},'description')
         if iscell(C{i}.description)
            fprintf(fileID,'# %s\n',C{i}.description{:}); % write multi-line description
         else
            fprintf(fileID,'# %s\n',C{i}.description); % write description
         end
         C{i}=rmfield(C{i},'description'); % do not write as field
      end % }}}
      % write PARM header {{{
      fprintf(fileID,' &%s\n',C{i}.header); % write header
      C{i}=rmfield(C{i},'header'); % do not write as field
      % }}}
      % write parameter fields and end section {{{
      if isfield(C{i},'N') % if diagnostic fields
         writediagfields(fileID,C{i}.N);
         C{i}=rmfield(C{i},'N'); % done with diag fields
      end
      writefields(fileID,C{i}); % write non-diagnostic fields
      fprintf(fileID,' &\n\n'); % end section
      % }}}
   end
   fclose(fileID);
end % }}}
function writefields(fileID,P) % {{{
   % WRITEFIELDS writes each field in P to fileID as a new line
   %  Each line takes the form ' fieldname=value,'
   fields=fieldnames(P);
   for i=1:length(fields)
      val=getfield(P,fields{i});
      fprintf(fileID,'  %s=%s,\n',fields{i},num2str(val));
   end
end  % }}}
function writediagfields(fileID,N) % {{{
   % WRITEDIAGFIELDS writes diagnostic fields for each output stream
   % INPUT: N is a struct. array with an element for each output stream n
   %  Each line takes the form fieldname(n)=value except for
   %  'fields' with take the form fields(1:length,n)='fields{1}', 'fields{2}', ...

   % loop over each output stream
   for n=1:length(N)
      subfields=fieldnames(N(n));
      % loop over each subfield
      for i=1:length(subfields)
         switch subfields{i}
            case 'fields'
               LHS=[subfields{i} '(1:' num2str(length(N(n).fields)) ',' num2str(n) ')'];
               dfields=getfield(N(n),subfields{i});
               dfields=strcat('''',dfields,''', ');
               RHS=strcat(dfields{:});
               fprintf(fileID,'  %s=%s\n',LHS,RHS); % write to file
            case 'levels'
               error('diag. levels not supported');
            otherwise
               LHS=[subfields{i} '(' num2str(n) ')'];
               RHS=num2str(getfield(N(n),subfields{i}));
               fprintf(fileID,'  %s=%s,\n',LHS,RHS); % write to file
         end
      end
      fprintf(fileID,'\n'); % line break
   end
end % }}}
