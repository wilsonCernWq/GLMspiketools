function [logli, dL, H] = Loss_GLM_logli_offset(prs,Xstruct)
% [logli, dL, H] = Loss_GLM_logli(prs,Xstruct)
%
% Compute negative log-likelihood of data undr the GLM model
% (with standard linear parametrization of stimulus kernel);
%
% Uses arbitrary nonlinearity 'nlfun' instead of exponential
%
% Inputs:
%   prs = [kprs - weights for stimulus kernel
%          dc   - dc current injection
%          ihprs - weights on post-spike current]
%
% Outputs:
%      logli = negative log likelihood of spike train
%      dL = gradient with respect to prs
%      H = hessian

% Unpack GLM prs;
nktot = Xstruct.nkx*Xstruct.nkt;   % total # params for k
kprs  = prs(1:nktot);
dc    = prs(nktot+1);
ihprs = prs(nktot+2:end-1);
bias  = prs(end);

% Extract some other stuff we'll use a lot
% matrix for interpolating from stimulus bins to spike train bins
M = Xstruct.Minterp;
% stimulus design matrix
XXstim = Xstruct.Xstim; 
% spike history design matrix
XXspik = Xstruct.Xsp;
% spikes
osps = Xstruct.osps; % original data
bsps = Xstruct.bsps; % binary spike vector
nsps = (sum(bsps));  % number of spikes
% number of bins in spike train vector
rlen = Xstruct.rlen;

% flag
ihflag = Xstruct.ihflag;
% absolute bin size for spike train (in sec)
dt = Xstruct.dtSp;
expO = exp(bias); % exp offset bias

% -------- Compute sum of filter reponses ---------------------------
if Xstruct.ihflag
   % stim-dependent + spikehist-dependent inputs
   
   Itot = M * (XXstim*kprs) + XXspik*ihprs + dc;
else
   % stim-dependent input only
   Itot = M * (XXstim*kprs) + dc;
end

% ---------  Compute output of nonlinearity  ------------------------
[f,df,ddf] = Xstruct.nlfun(Itot);

% ---------  Compute log-likelihood ---------------------------------
Trm1 =  sum(f + expO) * dt;  % non-spike term
Trm2 = -sum(log(f(bsps) + expO)); % spike term
logli = Trm1 + Trm2;

% ---------  Compute Gradient -----------------
if (nargout > 1)
   
   % Non-spiking terms (Term 1)
   dLdk0 = (df' * M * XXstim)';
   dLdb0 = sum(df);
   if ihflag
      dLdh0 = (df' * XXspik)';
   end
   dLdO0 = rlen * expO;
   
   % Spiking terms (Term 2)
   Msp = M(bsps,:); % interpolation matrix just for spike bins
   Nsp = osps(bsps);
   
   frac1 = Nsp .* df(bsps) ./ (f(bsps) + expO);
   dLdk1 = ((frac1' * Msp) * XXstim)';
   dLdb1 = sum(frac1);
   if ihflag
      dLdh1 = (frac1' * XXspik(bsps,:))';
   end
   dLdO1 = sum((Nsp .* expO) ./ (f(bsps) + expO));
   
   % Combine Term 1 and Term 2
   dLdk = dLdk0 * dt - dLdk1;
   dLdb = dLdb0 * dt - dLdb1;   
   if ihflag
      dLdh = dLdh0 * dt - dLdh1;
   else
      dLdh = [];
   end
   dLdO = dLdO0 * dt - dLdO1;
   dL = [dLdk; dLdb; dLdh; dLdO];
   
end

% ---------  Compute Hessian -----------------
if nargout > 2
   
   % --- Non-spiking terms -----
   % multiply each row of M with drr
   ddrrdiag = spdiags(ddf,0,rlen,rlen);
   % this is MUCH faster than using bsxfun, due to sparsity!
   ddqqIntrp = ddrrdiag*M; 
   
   % k and b terms
   Hk = (XXstim'*(M'*ddqqIntrp)*XXstim)*dt; % Hkk (k filter)
   Hb = sum(ddf)*dt;                      % Hbb (constant b)
   Hkb = (sum(ddqqIntrp,1)*XXstim)'*dt;   % Hkb (cross-term)
   if ihflag  % h terms
      Hh = (XXspik'*(bsxfun(@times,XXspik,ddf)))*dt;	% Hh (h filter)
      % (here bsxfun is faster than diagonal multiplication)
      
      Hkh = ((XXspik'*ddqqIntrp)*XXstim*dt)';	% Hhk (cross-term)
      Hhb = (ddf'*XXspik)'*dt;                  % Hhb (cross-term)
   else
      Hkh=[]; Hhb=[]; Hh=[];
   end
   
   % --- Add in spiking terms ----
   % needed weights from derivation of Hessian
   frac2 = (f(bsps).*ddf(bsps) - df(bsps).^2)./f(bsps).^2; 
   % rows of Msp re-weighted by these weights
   fr2Interp = spdiags(frac2,0,nsps,nsps)*Msp; 
   
   % Spiking terms, k and b
   Hk= Hk - XXstim'*(Msp'*fr2Interp)*XXstim; % Hkk (k filter)
   Hb =  Hb-sum(frac2);           % Hbb (constant b)
   Hb1 = sum(fr2Interp,1)*XXstim;  % Spiking term, k and const
   Hkb = Hkb - Hb1';
   if Xstruct.ihflag
      XXrrsp = bsxfun(@times,XXspik(bsps,:),frac2);
      Hh= Hh - XXspik(bsps,:)'*XXrrsp;
      % Const by h
      Hhb1 = sum(XXrrsp,1)';
      Hhb = Hhb - Hhb1;
      % k by h term
      Hkh0 = XXspik(bsps,:)'*(fr2Interp*XXstim);
      Hkh = Hkh-Hkh0';
   end
      
   HOk = (M * XXstim)' * (-osps .* df .* expO ./ ((f + expO).^2));
   HOb = (-osps .* df .* expO ./ ((f + expO).^2))' * df;
   if Xstruct.ihflag
      HOh = XXspik' * (-osps .* df .* expO ./ ((f + expO).^2));
   else
      HOh = [];
   end
   HOO = sum(-osps .* (expO ./ (f + expO) - expO * expO ./ ((f + expO).^2)) - expO * dt);

   HOcol = [HOk; HOb; HOh; HOO];
   HOrow = [HOk; HOb; HOh]';
   
   H = [[Hk Hkb Hkh]; [Hkb' Hb Hhb']; [Hkh' Hhb Hh]; HOrow];
   H = [H HOcol];
   
end

