function [prs0,Xstruct] = GLMsetupFit(gg, Stim)
%  [prs0,Xstruct] = setupfitting_GLM(gg, Stim)
%
%  Set initial parameters and build design matrix structure for GLM fitting
%
%  Inputs:
%     gg = glm param structure
%     Stim = stimulus (time along columns, other dims along rows)
%     maxsize = maximum # floats to store in design matrix
%                     (Affects # of chunks used to compute likelihood)
%  Output:
%     prs0 = initial parameters extracted from gg
%     Xstruct = struct with design matrix for stim and spike-history terms


% Initialize optimization param structure

% ---- Create struct and make stimulus design matrix ---------------------
if strcmp(gg.ktype, 'linear') % standard GLM
   
   % create design matrix structure
   Xstruct = GLMgg2xstruct(gg,Stim);  
   
   % extract parameter vector
   prs0 = [gg.kt(:); gg.dc; gg.ihw(:); gg.ihw2(:)];
   
elseif strcmp(gg.ktype, 'bilinear') % bilinearly-parametrized stim filter
   
   % create design matrix structure
   Xstruct = GLMgg2xstruct_bi(gg,Stim); 
   
   % extract parameter vector
   prs0 = [gg.kt(:); gg.kx(:); gg.dc; gg.ihw(:); gg.ihw2(:)];
   
elseif strcmp(gg.ktype, 'nobasis') 
      
   % create design matrix structure
   Xstruct = GLMgg2xstruct_nobasis(gg,Stim);
   
   % extract parameter vector
   prs0 = [gg.k(:); gg.dc; gg.ihw(:); gg.ihw2(:)];
   
elseif strcmp(gg.ktype, 'offset')
   
   % create design matrix structure
   Xstruct = GLMgg2xstruct_offset(gg,Stim);
   
   % extract parameter vector
   prs0 = [gg.kt(:); gg.dc; gg.ihw(:); gg.ihw2(:); gg.offset];
   
else
   
   error('unknown filter type (allowed types are ''linear'' or ''bilinear'')');
   
end

% ---- Make spike-history design matrix -----------------------------------
Xstruct = GLMinitSpikeHistDesignMat(gg,Xstruct);

% set nonlinearity
Xstruct.nlfun = gg.nlfun;

% compute mask (time bins to use for likelihood)
Xstruct.bmask = GLMinitFitMask(gg.mask,Xstruct.dtSp,Xstruct.rlen);

end

