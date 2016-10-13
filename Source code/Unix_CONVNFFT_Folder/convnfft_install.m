function convnfft_install
% function convnfft_install
% Installation by building the C-mex file needed for convnfft
%
% Author: Qi Dai <daiailiu04@yahoo.com>
% History
%  Original: 11/May/2015

arch=computer('arch');
mexopts = {'-O' '-v' ['-' arch]};
% 64-bit platform
if ~isempty(strfind(computer(),'64'))
    mexopts(end+1) = {'-largeArrayDims'};
end

% invoke MEX compilation tool
mex(mexopts{:},'inplaceprod.c');