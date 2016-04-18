function Xstruct = GLMgg2xstruct_nobasis(gg,Stim)
% Xstruct = initfit_stimDesignMat(gg,Stim)
%

% Initialize parameters relating to stimulus design matrix
%
% ---- Set up filter and stim processing params -------------------
% @remark
%   the initial value of K is still guessed using basis functions
% * size of kernel parameters
nkx = size(gg.k,2);
nkt = size(gg.k,1);
ncols = nkx * nkt; % total number of columns in stim design matrix
% * size of Stim matrix
[slen,swid] = size(Stim); % size of stimulus
upsampfactor = (gg.dtStim/gg.dtSp); % number of times by which spikes more finely sampled than stim
rlen = slen * upsampfactor; % length of spike train vector
%
% ---- Check size of filter and width of stimulus ----------
assert(nkx == swid,'Mismatch between stim width and kernel width');
%
% ---- Convolve stimulus with spatial and temporal bases -----
Xstruct.Xstim = zeros(slen,ncols);
basis = fliplr(eye(nkt));
for i = 1:nkx
   for j = 1:nkt
      Xstruct.Xstim(:,(i-1) * nkt + j) = ...
         sameconv(Stim(:,i),basis(:,j));     
   end
end
%
% ---- Set fields of Xstruct -------------------------------------
% nkx             : number stimulus spatial stimulus pixels
% nkt             : number of time bins in stim filter
% slen            : Total stimulus length (coarse bins)
% rlen            : Total spike-train bins (fine bins)
% upsamplefactor  : rlen / slen
% Minterp         : rlen x slen matrix for upsampling and downsampling
% dtStim          : time bin size for stimulus
% dtSp            : time bins size for spike train
Xstruct.nkx = nkx;
Xstruct.nkt = nkt;
Xstruct.slen = slen;
Xstruct.rlen = rlen;
Xstruct.upsampfactor = upsampfactor;
Xstruct.Minterp      = kron(speye(slen),ones(upsampfactor,1));
Xstruct.dtStim       = gg.dtStim;
Xstruct.dtSp         = gg.dtSp;
end