function [f,PSD] = readPSDfile(fname_fmt,ind,varargin)
%% function to read PSD file

if numel(varargin) < 1 %not sepcified
    avg_flag=false;
    zonal_flag=false;
elseif numel(varargin) == 1 
    avg_flag=varargin{1};
    zonal_flag = false;
else
    avg_flag=varargin{1};
    zonal_flag =varargin{2};
end


if avg_flag
    fname = sprintf(fname_fmt,ind);
    if zonal_flag
        A = readmatrix(fname,'NumHeaderLines',5);
    else
        A = readmatrix(fname,'NumHeaderLines',4);
    end

    f = A(:,1);
    PSD1 = A(:,2);
    fname = sprintf(fname_fmt,ind+111);
    if zonal_flag
        A = readmatrix(fname,'NumHeaderLines',5);
    else
        A = readmatrix(fname,'NumHeaderLines',4);
    end
    PSD2=A(:,2);
    PSD = (PSD1+PSD2)./2;

else
    fname = sprintf(fname_fmt,ind);
    if zonal_flag
        A = readmatrix(fname,'NumHeaderLines',5);
    else
        A = readmatrix(fname,'NumHeaderLines',4);
    end
    f = A(:,1);
    PSD = A(:,2);
end
end