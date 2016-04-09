function [STA,STC,RawMu,RawCov] = ...
   simpleSTC(stim,spik,nkt,maxsize,spikType)
% [STA,STC,RawMu,RawCov] = simpleSTC(Stim,sp,nkt,maxsize,spikType)
%
% Computes average and covariance of spike-triggered (or response weighted)
% stimuli and raw stimuli
%
% INPUT:
%    stim [N x M]    - stimulus matrix
%                       >> 1st dimension = time
%                       >> 2nd dimension = space
%    spik [N x 1]    - column vector of spike count in each time bin, which
%                       can be sparse
%    (or) [Nsp x 1]  - vector of spike times
%     nkt [1 x 1]    - No. of time samples to consider to be part of the 
%                       stimulus
% MaxSize [1 x 1] (optional)  - max No. of floats to store while computing
%                                cov (smaller = slower, but less memory 
%                                required)
% spikType string (optional)  - (added by WU,Qi) spercify types of spik 
%                                vector. 'spikeTimes' force the program
%                                to convert spik vector from time-series
%                                format into spike count format
%
%  OUTPUT:
%    STA [nkt x M]         - spike-triggered average (reshaped as matrix)
%    STC [nkt*M x nkt*M]   - spike-triggered covariance (covariance around 
%                             the mean)
%  RawMu [nkt*M x 1]       - mean of raw stimulus ensemble
% RawCov [nkt*M x nkt*M]   - covariance of raw stimulus ensemble
%
%  Notes:
%  (1) Ignores response before nkt bins
%  (2) Faster if only 2 output arguments (raw mean and covariance not 
%     computed)
%  (3) Reduce 'maxsize' if getting "out of memory" errors
%  (4) If inputs are not spike times or counts, response-weighted 
%     covariance may be more appropriate. See  "simpleRWC.m".
%
%  --------
%  Details:
%  --------
%   Let X = "valid" design matrix (from makeStimRows)
%       Y = sp (spike train)
%     nsp = sum(Y); % number of spikes (sum of responses)
%       N = size(Stim,1); % number of stimuli
%
%   then  STA = X'*Y / nsp
%         STC = X'*(X.*repmat(Y,1,ncols)))/(nsp-1) - STA*STA'*nsp/(nsp-1)
%         RawMu = sum(X)/N;
%         RawCov = X*X'/(N-1);
%
% Copyright 2010 Pillow Lab. All rights reserved.
% $Id$

%-------- Parse inputs  ---------------------------
% max chunk size; decrease if getting "out of memory"
if ~exist('maxsize','var')
   maxsize = 1e9;
end

[slen,swid] = size(stim); % stimulus size (# time bins x # spatial bins).
nsp = length(spik);

% Convert list of spike times to spike-count vector, if necessary
if ~exist('spikType', 'var') || strcmp(spikType,'spikeTimes')
   if (nsp ~= slen)
      fprintf(1, 'simpleSTC: converting spike times to counts\n');
      spik = hist(spik,1:slen)';
   end
end

spik(1:nkt-1) = 0;  % Remove spikes before time n
spsum = sum(spik);  % Sum of spikes
if (spsum == 0)
   spsum = 2;
end
nspnds = sum(spik~=0); % number of non-zero spike count

if nargout <= 2
   % ---------------------------------------------------
   % 1. Compute only the spike-triggered STA and STC

   Msz = nspnds * swid * nkt; % size of design matrix for spiking stimuli
   iisp = find(spik); % indices of spikes
   splen = length(iisp);
   
   if Msz < maxsize  % Compute in one chunk if small enough
      
      SS = makeStimRows(stim,nkt,iisp);
      STA = (SS'*spik(iisp))/spsum;
      if nargout > 1
         STC = SS'*bsxfun(@times,SS,spik(iisp))/(spsum-1)...
            - STA*STA'*spsum/(spsum-1);
      end
      
   else % Compute in multiple chunks if too large
      
      nchunk = ceil(Msz/maxsize);
      chunksize = ceil(length(iisp)/nchunk);
      fprintf(1,'simpleSTC: using %d chunks to compute STA/STC\n',nchunk);
      
      % Initialize on 1st chunk
      i0 = 1;
      imx = chunksize;
      SS = makeStimRows(stim,nkt,iisp(i0:imx));
      STA = (spik(iisp(i0:imx))'*SS)';
      if nargout > 1
         STC = SS'*bsxfun(@times,SS,spik(iisp(i0:imx)));
      end
      
      % Compute for remaining chunks
      for j = 2:nchunk
         i0 = chunksize*(j-1)+1;
         imx = min(chunksize*j, splen);
         SS = makeStimRows(stim,nkt,iisp(i0:imx));
         STA = STA + (spik(iisp(i0:imx))'*SS)';
         if nargout > 1
            STC = STC + SS'*bsxfun(@times,SS,spik(iisp(i0:imx)));
         end
      end
      
      % normalize by number of samples
      STA = STA/spsum;
      if nargout > 1
         STC = STC/(spsum-1) - STA*STA'*spsum/(spsum-1);
      end
      
   end
   
else
   % ---------------------------------------------------
   % 2. Compute both the spike-triggered and raw stimulus means and 
   % covariances

   spik = spik(nkt:end);    % Ignore spikes before time bin nkt
   slen = length(spik);   % length of time indices for stim and spikes
   Msz = slen*swid*nkt;   % Size of full stimulus matrix
   
   if Msz < maxsize  % Check if stimulus is small enough to do in one chunk
      
      SS = makeStimRows(stim,nkt,'valid'); % Convert stimulus to matrix where each row is one stim
      
      % Compute raw mean and covariance
      RawMu = mean(SS)';
      RawCov = (SS'*SS)/(slen-1)-RawMu*RawMu'*slen/(slen-1);
      
      % Compute spike-triggered mean and covariance
      iisp = find(spik>0);
      spvec = spik(iisp);
      STA = (spvec'*SS(iisp,:))'/spsum;
      STC = SS(iisp,:)'*bsxfun(@times,SS(iisp,:),spvec)/(spsum-1) ...
         - STA*(STA'*spsum/(spsum-1));
      
   else  % Compute Full Stim matrix in chunks, compute mean and cov on chunks
      
      nchunk = ceil(Msz/maxsize);
      chunksize = ceil(slen/nchunk);
      fprintf(1, 'simpleSTC: using %d chunks to compute covariance\n', nchunk);
      
      % Compute statistics on first chunk
      SS = makeStimRows(stim(1:chunksize+nkt-1,:),nkt,'valid');  % convert stimulus to "full" version
      spvec = spik(1:chunksize);
      iisp = find(spvec>0);
      RawMu = sum(SS)';
      RawCov = SS'*SS;
      STA = (spvec(iisp)'*SS(iisp,:))';
      STC = SS(iisp,:)'*bsxfun(@times,SS(iisp,:),spvec(iisp));
      
      % add to mean and covariance for remaining chunks
      for j = 2:nchunk;
         i0 = chunksize*(j-1)+1;  % starting index for chunk
         imax = min(slen,chunksize*j);  % ending index for chunk
         SS = makeStimRows(stim(i0:imax+nkt-1,:),nkt,'valid');
         spvec = spik(i0:imax);
         iisp = find(spvec);
         
         RawMu = RawMu + sum(SS)';
         RawCov = RawCov +SS'*SS;
         STA = STA + (spvec(iisp)'*SS(iisp,:))';
         STC = STC + SS(iisp,:)'*bsxfun(@times,SS(iisp,:),spvec(iisp));
      end
      
      % divide means and covariances by number of samples
      RawMu = RawMu/slen;
      RawCov = RawCov/(slen-1) - RawMu*RawMu'*slen/(slen-1);
      
      STA = STA/spsum;
      STC = STC/(spsum-1) - STA*STA'*spsum/(spsum-1);
      
   end
   
end

STA = reshape(STA,[],swid);  % reshape to have same width as stimulus
