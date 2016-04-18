function gg = GLMxstruct2gg(gg,prs,Xstruct)
% gg = reinsertFitPrs_GLM(gg,prs,Xstruct)
%
% DESCRIPTION
%
% After fitting, reinsert params into GLM param structure
%
if strcmp(gg.ktype, 'linear')
   
   % Extract relevant size information
   nkt = Xstruct.nkt;
   nkx = Xstruct.nkx;
   nktot = nkt * nkx;
   nh = Xstruct.nh;
   nh2 = Xstruct.nh2;
   
   % Insert params
   gg.kt = reshape(prs(1:nktot),nkt,nkx);
   gg.k = gg.ktbas * gg.kt;
   gg.dc = prs(nktot+1);
   gg.ihw = reshape(prs(nktot+2:nktot+nh+1), nh,1);
   gg.ihw2 = reshape(prs(nktot+nh+2:end), nh2, []);
   
elseif strcmp(gg.ktype, 'nobasis')
   
   % Extract relevant size information   
   nkt = Xstruct.nkt;  
   nkx = Xstruct.nkx;
   nktot = nkt*nkx;
   nh = Xstruct.nh;  
   nh2 = Xstruct.nh2;

   % Insert params
   basis = fliplr(eye(nkt));   
   gg.k  = basis * reshape(prs(1:nktot),nkt,nkx);
   gg.dc = prs(nktot+1);
   gg.ihw = reshape(prs(nktot+2:nktot+nh+1), nh,1);
   gg.ihw2 = reshape(prs(nktot+nh+2:end), nh2, []);

elseif strcmp(gg.ktype, 'bilinear')
   
   % Put returned vals back into param structure ------
   krank = Xstruct.krank;
   nktprs = Xstruct.nkt*krank;
   nkxprs = Xstruct.nkx*krank;
   nktot = nktprs+nkxprs;
   nh = Xstruct.nh;
   nh2 = Xstruct.nh2;
   
   % Insert params into struct
   gg.kt = reshape(prs(1:nktprs),[],krank);
   gg.kx = reshape(prs(nktprs+1:nktprs+nkxprs),[],krank);
   gg.dc = prs(nktot+1);
   gg.ihw = reshape(prs(nktot+2:nktot+nh+1), nh,1);
   gg.ihw2 = reshape(prs(nktot+nh+2:end), nh2, []);
   gg.k = (gg.ktbas*gg.kt)*(gg.kxbas*gg.kx)';
   
end

% Ensure these bases are 'empty' if no spike history filters
if isempty(gg.ihw)
   gg.ihbas = [];
end
if isempty(gg.ihw2)
   gg.ihbas2 = [];
end

% Insert spike-history filters
gg.ih = [gg.ihbas*gg.ihw, gg.ihbas2*gg.ihw2];
end
