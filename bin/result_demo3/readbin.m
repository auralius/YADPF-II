function output = readbin(bin_file, n_channels, data_type)
% This function reads a general binary file that contains a sequence of
% numbers. This function then reshapes those numbers into a matrix of
% N x n_channels where N depends on how many number are there.
%
% This function has been tested for binary files created by LabVIEW and
% Armadillo library.
%
% bin_file: the binary file must be composed by n_channels of double or
% int32 data type.
%
% n_channels: number of channels
% output: channel-1 is in row-1
%         channel-2 is in row-2
%         channel-3 is in row-3
%         and so on ...
%
% Example = data = readbin('logged_data.bin', 45, 'double');
%
% Contact: manurunga@yandex.com
    
    if nargin < 3
        data_type = 'double';
    end

    fid = fopen(bin_file, 'r');

    % Find out the length of the file in bytes
    fseek(fid, 0, 'eof');
    filelength = ftell(fid);
    frewind(fid);
    
    % Read the file:
    % 1 double = 8 bytes, read per 8 bytes, convert into double
    % 1 int32 = 4 bytes, read per 4 bytes, convert into int32
    if strcmp(data_type, 'int32')
        output = fread(fid, filelength/4, 'int32');
    elseif strcmp(data_type, 'double')
        output = fread(fid, filelength/8, 'double');
    elseif strcmp(data_type, 'int64')
        output = fread(fid, filelength/8, 'int64');
    else
        error('Define correct data type: int32 or double?')
    end
    
    % Length of the output must be factor of n_channels
    % Reshape accordaingly to the number of n_channels
    if mod(length(output)/n_channels, 1) == 0 % Integer value?
        output = reshape(output, n_channels, length(output)/n_channels);
    else
        disp('Are you sure you gave the right number of channels?')
        output = 0;
    end
    
    fclose(fid);
end