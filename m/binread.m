function D=binread(fname,prec,arrsize)
% BINREAD read data from binary file into a matlab array D.
% Assumes big-endian architecture, and given precision
% and array size.
%
% fname: filename or path (string)
% prec: 4 or 8 for number of bits
% arrsize: dimensions of D (array)
%
% D: array of requested dimension
%
% Example:
% D=binread('bathy.bin',8,[Nx,Ny])

fid=fopen(fname,'r','b');
switch prec
	case 8
		D=fread(fid,inf,'real*8');
	case 4
		D=fread(fid,inf,'real*4');
	otherwise
		error('give precision of data');
end
D=reshape(D,arrsize);
fclose(fid);
