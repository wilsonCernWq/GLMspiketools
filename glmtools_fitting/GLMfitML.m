function [gg,fval,H,Xstruct] = GLMfitML(gg,Stim,optimArgs)
%  [gg,fval,H,Xstruct] = MLfit_GLM(gg,Stim,optimArgs)
%
%  Computes the ML estimate for GLM params, using grad and hessians.
%  Assumes basis for temporal dimensions of stim filter
%
%  Inputs:
%     gg = param struct
%     Stim = stimulus
%     optimArgs = cell array of optimization params (optional)
%
%  Outputs:
%     ggnew = new param struct (with ML params);
%     fval = negative log-likelihood at ML estimate
%        H = Hessian of negative log-likelihood at ML estimate
%  Xstruct = structure with design matrices for spike-hist and stim terms

% Set optimization parameters
if nargin > 2
   opts = optimset('Gradobj','on','Hessian','on', optimArgs{:});
else
   opts = optimset('Gradobj','on','Hessian','on','display','iter');
end

% --- Create design matrix extract initial params from gg ----------------
[prs0,Xstruct] = GLMsetupFit(gg,Stim);

% --- Set loss function --------------------------------------------------
if strcmp(gg.ktype, 'linear') || strcmp(gg.ktype, 'nobasis')
   
   if isequal(Xstruct.nlfun,@expfun) || isequal(Xstruct.nlfun,@exp)
      % loss function for exponential nonlinearity
      lfunc = @(prs)Loss_GLM_logli_exp(prs,Xstruct); 
   else
      % loss function for all other nonlinearities
      lfunc = @(prs)Loss_GLM_logli(prs,Xstruct); 
   end
   
   % --- minimize negative log likelihood --------------------------------
   [prsML,fval] = fminunc(lfunc,prs0,opts); % find ML estimate of params
   
   % Compute Hessian if desired
   if nargout > 2
      [fval,~,H] = Loss_GLM_logli(prsML,Xstruct);
   end

elseif strcmp(gg.ktype, 'bilinear');      
   
   % loss function for exponential nonlinearity
   lfunc = @(prs)Loss_GLM_logli_bi(prs,Xstruct);
   
   % optimize negative log-likelihood for prs
   [prsML,fval] = fminunc(lfunc,prs0,opts); 
   
   % Compute Hessian at maximum, if requested
   if nargout > 2
      [fval,~,H] = Loss_GLM_logli_bi(prsML);
   end
   
end

% Put returned vals back into param structure
gg = GLMxstruct2gg(gg,prsML,Xstruct);

% ----------------------------------------------------
% Optional debugging code
% ----------------------------------------------------
%
% ------ Check analytic gradients and Hessians -------
%  HessCheck(lfunc,prs0,opts);
%  HessCheck_Elts(@Loss_GLM_logli, [1 12],prs0,opts);
%  tic; [lival,J,H]=lfunc(prs0); toc;

