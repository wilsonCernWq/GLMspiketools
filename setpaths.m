% SETPATHS.m - GLMspiketools code repository
%
% This simple script sets the path to include relevant directories for the
% GLMspiketools code package.  You must 'cd' into this directory in order
% to evaluate it.
%
% More info: http://pillowlab.princeton.edu/code_GLM.html
% Github page: https://github.com/pillowlab/GLMspiketools

basedir = pwd;  % The directory where this script lives
GLMspiketools_path = ...
   {[basedir '/glmtools_fitting/'],...
    [basedir '/glmtools_fitting/full'],...
    [basedir '/glmtools_fitting/bilinear'],...
    [basedir '/glmtools_fitting/nobasis'],...
    [basedir '/glmtools_fitting/offset'],...
    [basedir '/glmtools_misc/'],...
    [basedir '/nlfuns/']};

% Add a bunch sub-directories (with absoluate path names)
for m_i = 1:length(GLMspiketools_path)
   addpath(GLMspiketools_path{m_i});
end
clear m_i;
