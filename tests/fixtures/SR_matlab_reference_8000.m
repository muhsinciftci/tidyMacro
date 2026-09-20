clear; close all; clc; warning off all
root = '/Users/muhsinciftci/Desktop/Packages/tidyMacro/inst/Bianchi';
out  = '/private/tmp/claude-501/-Users-muhsinciftci-Desktop-Packages-tidyMacro/1f6bfb37-4fd7-4905-bae3-0f37fc977745/scratchpad';
addpath(fullfile(root,'VAR')); addpath(fullfile(root,'Auxiliary'));
addpath(fullfile(root,'Stats')); addpath(fullfile(root,'Utils'));

ND = 8000;
rng(42,'twister');
raw = readcell(fullfile(root,'Replic','Uhlig2005','Uhlig2005_Data.xlsx'),'Sheet','Sheet1');
mnem = raw(2,2:end); data = cellfun(@double, raw(3:end,2:end));
for ii=1:numel(mnem); DATA.(mnem{ii}) = data(:,ii); end
for v = {'y','pi','comm','nbres','res'}; DATA.(v{1}) = 100*DATA.(v{1}); end
X = nan(size(data,1),numel(mnem));
for ii=1:numel(mnem); X(:,ii) = DATA.(mnem{ii}); end
SIGN = [0 0 0 0 0 0; -1 0 0 0 0 0; -1 0 0 0 0 0; 0 0 0 0 0 0; -1 0 0 0 0 0; 1 0 0 0 0 0];
VARopt = VARoption; VARopt.mnem = mnem;
VARopt.nsteps=60; VARopt.ndraws=ND; VARopt.inference=1; VARopt.sr_hor=6;
VARopt.pctg=68; VARopt.mult=1000; VARopt.ident='sign'; VARopt.R=SIGN;
S = VARmodel(X,12,1,VARopt);
writematrix([S.IRmed(:,:,1) S.IRinf(:,:,1) S.IRsup(:,:,1)], fullfile(out,'ml_uhlig_irf_hi.csv'));
fprintf('UHLIG HI done accept=%.6f\n', S.accept_rate);
fprintf('ALL DONE HI\n');
