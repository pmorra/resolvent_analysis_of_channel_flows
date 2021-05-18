function [PHI,MU] = eigW(Ptot,W)

[PHI,LMU] = eig(Ptot*W); MU = abs(diag(LMU));
[MU,ii] = sort(MU,'descend');
PHI = PHI(:,ii);
rk = rank(Ptot*W);
MU = MU(1:rk); 
for isp = 1:rk
  thenorm = sqrt(PHI(:,isp)'*W*PHI(:,isp));
  PHI(:,isp) = PHI(:,isp)/thenorm;
end
PHI = PHI(:,1:rk);
end