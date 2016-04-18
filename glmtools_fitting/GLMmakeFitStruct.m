function gg = GLMmakeFitStruct(varargin)
% gg = makeFittingStruct_GLM
%      (fittype,dtStim,dtSp,klength,nkbasis,k0,nhbasis,lasthpeak)
%
% DESCRIPTION
%
%  Initialize parameter structure for fitting GLM, with normal 
%  parametrization of stim kernel
%
% Inputs:
%
%     dtStim = bin size of stimulus (in s)
%     dtSpik = bin size for spike train (in s)
%    klength = temporal length of stimulus filter, in # of bins (optional)
%    nkbasis = # temporal basis vectors for stim filter         (optional)
%    nhbasis = # temporal basis vectors for spike hist filter h (optional)
%  lasthpeak = time of peak of last h basis vector, in s        (optional)
%         k0 = initial estimate of stimulus filter              (optional)   
%
% Outputs:
%         gg = GLM-fitting param structure

% convert input arguments
fittype = varargin{1};
dtStim  = varargin{2};
dtSpik  = varargin{3};
if nargin > 3
   klength = varargin{4};
   nkbasis = varargin{5};
   if nargin > 5
      k0   = varargin{6};
   end   
   if nargin > 6
      hlength = varargin{8};   
      nhbasis = varargin{7};
   end   
end
% ---- Set up fitting structure -------------------------------
%
% basic
gg.k = [];  % Actual stim filter k
gg.ih = []; % Actual post-spike filter ih
gg.dc = 0;  % "dc" or constant input (determines baseline spike rate)
gg.nlfun = @expfun; % default nonlinearity: exponential
%
% post-spike history kernel
gg.ihw = [];  % h filter weights 
gg.iht = [];  % h filter time points
gg.ihbas = []; % basis for h filter
gg.ihbasprs = []; % parameters for basis for h-filter
gg.ihw2 = [];      % weights for coupling filters 
gg.ihbas2 = [];    % basis for coupling filters
gg.ihbasprs2 = []; % parameters for coupling filter basis 
%
% stimuli kernel
gg.kt = [];       % basis weights for stimulus filter k
gg.ktbas = [];    % temporal basis for stimulus filter k
gg.ktbasprs = []; % parameters for basis for k-filter
%
% spike vector
gg.sps = [];  % spike times (in s)
gg.sps2 = [];      % spike times of coupled neurons
%
% temporal bin size
gg.dtStim = dtStim;  % time bin size for stimulus 
gg.dtSp = dtSpik;    % time bin for spike train 
%
% others
gg.mask = []; % list of intervals to ignore when computing likelihood
gg.couplednums = []; % numbers of coupled cells

% ----- Set up temporal basis for stimulus kernel -----------
if nargin > 3
   
   assert((klength>nkbasis), ...
      'klength should be bigger than number of temporal basis vectors');
   
   % number of "identity" basis vectors
   ktbasprs.neye = 0;
   % Number of raised-cosine vectors to use
   ktbasprs.ncos = nkbasis;
   % Position of 1st and last bump
   ktbasprs.kpeaks = [0 klength*(1 - 1.5/nkbasis)];
   % Offset for nonlinear scaling (larger -> more linear)
   ktbasprs.b = 10;
   
   [~,ktbasis] = makeBasis_StimKernel(ktbasprs,klength);
   gg.ktbas = ktbasis;
   gg.ktbasprs = ktbasprs;
   
end

if (nargin > 5) && (~isempty(k0))
   
   % initialize k filter in this basis
   gg.kt = (gg.ktbas'*gg.ktbas)\(gg.ktbas'*k0);
   gg.k = gg.ktbas*gg.kt;
   
end

% ----- Set up basis for post-spike filter -----------------------
if nargin > 6
   
   % number of basis vectors for post-spike kernel
   ihbasprs.ncols = nhbasis;
   % peak location for first and last vectors
   ihbasprs.hpeaks = [dtSpik hlength];
   % how nonlinear to make spacings (rough heuristic)
   ihbasprs.b = hlength/5;
   % absolute refractory period (optional)
   ihbasprs.absref = [];
   
   [iht,ihbas] = makeBasis_PostSpike(ihbasprs,dtSpik);
   gg.iht = iht;
   gg.ihbas = ihbas;
   gg.ihbasprs = ihbasprs;
   gg.ihw = zeros(size(ihbas,2),1);
   % Set ih to be the current value of ihbas*ihw
   gg.ih = gg.ihbas*gg.ihw;
   
end

% specify type of parametrization of filter ('linear' vs. 'bilinear')
if ~exist('fittype','var') || strcmp(fittype,'linear')
   
   gg.ktype = 'linear';

elseif strcmp(fittype,'bilinear')
   
   % input argument
   krank = varargin{9};
   
   % Set additional fields needed by bilinear GLM
   gg.kx = [];
   gg.kxbas=[];
   gg.kbasprs = [];
   gg.krank = krank;
   gg.ktype = 'bilinear';
   
   % if initial filter passed in, use svd to set up initial K params 
   % (bilinearly parametrized)
   if (length(varargin) > 5) && (~isempty(varargin{3}))
      
      [u,s,v] = svd(k0);
      mux = mean(k0)';
      nkt = size(gg.ktbas,2);
      nkx = size(k0,2);
      kt = zeros(nkt,krank);
      kx = zeros(nkx,krank);
      for j = 1:krank;
         if v(:,j)'*mux < 0  % Flip sign if shape looks opposite to k0
            u(:,j) = -u(:,j); v(:,j) = -v(:,j);
         end
         kt(:,j) = (gg.ktbas'*gg.ktbas)\(gg.ktbas'*(u(:,j)*s(j,j)));
         kx(:,j) = v(:,j);
      end
      gg.kxbas = speye(nkx);
      gg.kt = kt;
      gg.kx = kx;
      gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';
      
   end
  
elseif strcmp(fittype,'nobasis')
   
   gg.ktype = 'nobasis';
   gg = rmfield(gg, {'ktbasprs','ktbas','kt'});
   
elseif strcmp(fittype,'offset')
   
   gg.ktype = 'offset';
   gg.offset = 0.0;
   
end

end

