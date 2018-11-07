function WindowedData = computeWinFunction(Data,WinSizeInSamps,functionHandle)
% Computes the windowed function for data X.
%
% @Input:
%   Data (required)
%     A nx1 numeric/logical column vector with the data. If this is a row
%     vector, a 'DJM:DimensionMismatch'-warning and the dimension is changed. If
%     this is a matrix, a 'DJM:DimensionMismatch'-error is thrown.
%
%   WinSizeInSamps (required)
%     An integer scalar >= 3 with the window size of the function. This has to
%     be odd (centered window), otherwise a 'DJM:WinSizeMustBeOdd'-warning will
%     be thrown and the window size will be increased by 1. Large window sizes
%     will cause significant memory load and may even crash the system.
%
%   functionHandle (required)
%     A functionhandle which evaluates column oriented data and generates a
%     scalar output (e.g. STD, MIN, ...).
%
% @Output:
%   WindowedData
%     The windowed data. To avoid start-/end-condition errors,
%     (WinSizeInSamps-1)/2 samples in the beginning and the end are set to NaN
%     (if the output of 'functionHandle' was numeric) or false (if it was
%     logical).
%
% @Remarks:
%   - No error checking is done on the type of 'Data' or 'functionHandle'.
%   - The function is not tested for complex inputs.
%   - If 'functionHandle' is @mean of @sum, convolution is used, which is again
%     faster then the reshaping approach. These shortcuts are inspired by Cris
%     Luengo's comment on Jeff Burkey's MOVINGSTATS at:
%       www.mathworks.com/matlabcentral/fileexchange/29029
%     However, the convolution results are only the same as the reshape approach
%     to a precision of roughly 10 digits. 
%
% @Dependancies
%   none
%
% @Changelog
%   2014-11-13 (DJM): Release.
%   2015-03-11 (DJM): Added logical vs. numeric support & function header.
%   2016-09-01 (DJM):
%     - Changed Id computation to use BSXFUN.
%     - Added CONV shortcuts for @sum & @mean which are again much faster.
%     - Added input argument checks for 'Data'.

%% Check input
%>Start: Comment this to run on any ML release.
if ~iscolumn(Data)
	if isrow(Data)
		warning('DJM:DimensionMismatch', ['Expected ''Data'' to be a column ' ...
			'vector. Got row instead!']);
		Data = Data';
	else % <=> ismatrix
		error('DJM:DimensionMismatch', '''Data'' must be a column vector!');
	end
end
%<end

if mod(WinSizeInSamps,2) == 1
	WinSizeHalfInSamps = floor(WinSizeInSamps/2);
else
	WinSizeHalfInSamps = WinSizeInSamps/2;
	WinSizeInSamps     = WinSizeInSamps+1;
	warning('DJM:WinSizeMustBeOdd', ['Window size has to be odd (%d). Using '...
		'%d instead!'], WinSizeInSamps-1, WinSizeInSamps);
end

%% Compute
nEnd = length(Data)-WinSizeInSamps+1;

%Compute the actual windowed function using...
if isequal(functionHandle,@mean)
	%...CONV shortcut for mean.
	WindowedData = conv(Data, repmat(1/WinSizeInSamps,WinSizeInSamps,1),'valid');
elseif isequal(functionHandle,@sum)
	%...CONV shortcut for sum.
	WindowedData = conv(Data, ones(WinSizeInSamps,1), 'valid');
else
	%...reshape approach using column-oriented matrix.
	WindowedData ...
		= functionHandle(Data(bsxfun(@plus,(1:WinSizeInSamps)',0:nEnd-1)))';
end

%Treat start-/end-condition based on type.
if islogical(WindowedData)
	invalid = @false;
else
	invalid = @NaN;
end
Dummy = invalid(WinSizeHalfInSamps,1);
WindowedData = [Dummy;WindowedData;Dummy];