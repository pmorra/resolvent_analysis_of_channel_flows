function [S,Var,psVar,f] = my_welch_for_Scov(F,nDFT,nOVL,dt,varargin)
%
% Computes the power spectral density matrix S = f*f', with
% S = S(x,x',omega) via the Welch Method.
%
% INPUT  F   : size(F) =  [Nt, Nx], with F the snapshots matrix
%        nDFT: number of points for the DFT on the bins (and number of
%              points per bin)
%        nOVL: number of overlapping points
%        dt  : sampling interval
%         
%        (Optional)
%        winExp : exponent for the windowing function sin(t)^winExp
%                 winExp = 2 is DEFAULT
%        diffWin: if true uses the first derivative of the window
%                 diffWin = false is DEFAULT
%        use_par: if true works activate parpool with the maxNumWorkers
%                 use_par = true is DEFAULT
%
% OUTPUT S     : size(S) = [nx,nx,nDFT], the spectral density matrix S
%        Var   : variance of each entry of S
%        psVar : pseudo-variance of each entry of S
%        f     : size(f) = [nDFT], vector of frequencies
%


% Set optional inputs
if nargin > 4; winExp = varargin{1}; else winExp = 2; end
if nargin > 5; diffWin = varargin{2}; else diffWin = false; end
if nargin > 6; use_par = varargin{3}; else use_par = true; end

% Extract dimensions
nx = size(F,2);
nt = size(F,1);
nBlks = floor((nt-nOVL)/(nDFT-nOVL));

% Obtain frequency axis
f = (0:nDFT-1)/dt/nDFT;
if mod(nDFT,2)==0
    f(nDFT/2+1:end) = f(nDFT/2+1:end)-1/dt;
else
    f((nDFT+1)/2+1:end) = f((nDFT+1)/2+1:end)-1/dt;
end
nFreq = length(f);

% Initialize containers for Welch
Q_hat = zeros(nFreq,nx,nBlks);
Q_blk       = zeros(nDFT,nx);

% Window for Welch
if ~diffWin
  % if the exponent is 2, this is Hann window
  dum = sin(pi*(0:nDFT-1)/nDFT).^winExp;
  window = dum.'; clear dum;  
  winWeight = 1/mean(window);
else
  dum = sin(pi*(0:nDFT-1)/nDFT).^winExp;
  winWeight = 1/mean(dum);    
  dum = (winExp)*sin(pi*(0:nDFT-1)/nDFT).^(winExp-1)...
        .*cos(pi*(0:nDFT-1)/nDFT)*pi/nDFT/dt;
  window = dum.'; clear dum; 
end

% Create sub-blocks
for iBlk = 1:nBlks
  % get time index for present block
  offset                  = min((iBlk-1)*(nDFT-nOVL)+nDFT,nt)-nDFT;
  timeIdx                 = (1:nDFT) + offset;
  % extract block
  Q_blk   = (F(timeIdx,:));
  Q_blk                   = bsxfun(@times,Q_blk,window);
  Q_blk_hat               = winWeight/nDFT*fft(Q_blk);
  Q_blk_hat               = Q_blk_hat(1:nFreq,:);
  Q_hat(:,:,iBlk)         = Q_blk_hat;
end

% Initialize containers for average over blocks
S = zeros(nx,nx,nFreq);
dVar = cell(nFreq,1);
dpsVar = cell(nFreq,1);

% Compute average over blocks
for iFreq = 1:nFreq
  Q_hat_f = squeeze(Q_hat(iFreq,:,:));
  if size(Q_hat_f,1) == 1; Q_hat_f = Q_hat_f.'; end
  S(:,:,iFreq) = Q_hat_f*Q_hat_f'/nBlks;
%   if ~use_par % [SERIAL] Compute variance and pseudo-variance over blocks
%     dVar{iFreq} = zeros(nx,nx);
%     dpsVar{iFreq} = zeros(nx,nx);
%     for iBlk = 1:nBlks
%       test = squeeze(squeeze(Q_hat(iFreq,:,iBlk))).';
%       dVar{iFreq} = dVar{iFreq} + abs(test*test').^2/nBlks;
%       dpsVar{iFreq} = dpsVar{iFreq} + (test*test').^2/nBlks;
%     end
%     dVar{iFreq} = dVar{iFreq} - abs(S(:,:,iFreq)).^2;
%     dpsVar{iFreq} = dpsVar{iFreq} - S(:,:,iFreq).^2;
%   end
end

% % [PARALLEL] Compute variance and pseudo-variance over blocks
% if use_par
%   mycl = parcluster('local');
%   parpool(mycl.NumWorkers);
%   parfor iFreq = 1:nFreq
%     dVar{iFreq} = zeros(nx,nx);
%     dpsVar{iFreq} = zeros(nx,nx);
%     for iBlk = 1:nBlks
%       test = squeeze(squeeze(Q_hat(iFreq,:,iBlk))).';
%       dVar{iFreq} = dVar{iFreq} + abs(test*test').^2/nBlks;
%       dpsVar{iFreq} = dpsVar{iFreq} + (test*test').^2/nBlks;
%     end
%     dVar{iFreq} = dVar{iFreq} - abs(S(:,:,iFreq)).^2;
%     dpsVar{iFreq} = dpsVar{iFreq} - S(:,:,iFreq).^2;
%   end
%   delete(gcp('nocreate'));
% end

Var = 0*S; psVar = 0*S;
% for iFreq = 1:nFreq
%   Var(:,:,iFreq) = dVar{iFreq};
%   psVar(:,:,iFreq) = dpsVar{iFreq};
% end

% Print out
fprintf(['\n', ...
'---- Welch function (Pierluigi) ----\n\n',...
'num Blocks = ',num2str(nBlks),'\n',...
'num DFT = ',num2str(nDFT),'\n',...
'\n------------------------------------\n',...
]);
end
