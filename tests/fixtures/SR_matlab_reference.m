% Ground-truth reference for the tidyMacro sign / narrative / sign+iv port.
% Runs VAR Toolbox 4.0 on the three replication datasets and writes numbers
% only.  Deliberately does NOT call SaveFigure or touch inst/Bianchi/Replic.

clear; close all; clc
warning off all

root = '/Users/muhsinciftci/Desktop/Packages/tidyMacro/inst/Bianchi';
out  = '/private/tmp/claude-501/-Users-muhsinciftci-Desktop-Packages-tidyMacro/1f6bfb37-4fd7-4905-bae3-0f37fc977745/scratchpad';
addpath(fullfile(root,'VAR'));
addpath(fullfile(root,'Auxiliary'));
addpath(fullfile(root,'Stats'));
addpath(fullfile(root,'Utils'));

% =====================================================================
% 1. UHLIG (2005) — sign restrictions
% =====================================================================
rng(42,'twister');
f   = fullfile(root,'Replic','Uhlig2005','Uhlig2005_Data.xlsx');
raw = readcell(f,'Sheet','Sheet1');
mnem = raw(2,2:end);
data = cellfun(@double, raw(3:end,2:end));
sc   = {'y','pi','comm','nbres','res'};
for ii = 1:numel(mnem)
    DATA.(mnem{ii}) = data(:,ii);
end
for ii = 1:numel(sc); DATA.(sc{ii}) = 100*DATA.(sc{ii}); end
X = nan(size(data,1), numel(mnem));
for ii = 1:numel(mnem); X(:,ii) = DATA.(mnem{ii}); end

SIGN = [ 0 0 0 0 0 0; -1 0 0 0 0 0; -1 0 0 0 0 0;
         0 0 0 0 0 0; -1 0 0 0 0 0;  1 0 0 0 0 0];

VARopt = VARoption;
VARopt.mnem = mnem;
VARopt.nsteps = 60; VARopt.ndraws = 500; VARopt.inference = 1;
VARopt.sr_hor = 6;  VARopt.pctg = 68;
VARopt.ident = 'sign'; VARopt.R = SIGN;
SRu = VARmodel(X, 12, 1, VARopt);

writematrix([SRu.IRmed(:,:,1) SRu.IRinf(:,:,1) SRu.IRsup(:,:,1)], ...
    fullfile(out,'ml_uhlig_irf.csv'));
writematrix(SRu.Bmed, fullfile(out,'ml_uhlig_Bmed.csv'));
writematrix([SRu.accept_rate SRu.ndraws_tried], fullfile(out,'ml_uhlig_diag.csv'));
fprintf('UHLIG done. accept_rate = %.6f\n', SRu.accept_rate);

% =====================================================================
% 2. ANTOLIN-DIAZ & RUBIO-RAMIREZ (2018) — sign vs sign + narrative
% =====================================================================
clearvars -except root out
rng(42,'twister');
f     = fullfile(root,'Replic','ADRR2018','ADRR2018_Data.xlsx');
raw   = readcell(f,'Sheet','Sheet1');
dates = raw(3:end,1);
mnem  = raw(2,2:end);
data  = cellfun(@double, raw(3:end,2:end));

SIGN = [ 0 0 0 0 0 0; -1 0 0 0 0 0; -1 0 0 0 0 0;
         0 0 0 0 0 0; -1 0 0 0 0 0;  1 0 0 0 0 0];

VARopt = VARoption;
VARopt.mnem = mnem;
VARopt.nsteps = 60; VARopt.ndraws = 500; VARopt.sr_hor = 6;
VARopt.pctg = 68;   VARopt.sr_draw = 500000; VARopt.inference = 1;
VARopt.mult = 100;  VARopt.dates = dates;

VARopt.ident = 'sign'; VARopt.R = SIGN;
SRa = VARmodel(data, 12, 0, VARopt);
fprintf('ADRR sign done. accept_rate = %.6f\n', SRa.accept_rate);

R.sign = SIGN;
R.narr_sign.shock = 1; R.narr_sign.period = '1979m10'; R.narr_sign.sign = 1;
R.narr_dom.shock  = 1; R.narr_dom.period  = '1979m10'; R.narr_dom.var   = 6;
VARopt.R = R;
NSRa = VARmodel(data, 12, 0, VARopt);
fprintf('ADRR narrative done. accept_rate = %.6f\n', NSRa.accept_rate);

% Normalisation of GO_ADRR2018.m: 25bp on ff at impact, applied to all draws
ffv = 6; mps = 1;
scS = (0.25 / median(squeeze(SRa.IRall(1,ffv,mps,:))))  * 100;
scN = (0.25 / median(squeeze(NSRa.IRall(1,ffv,mps,:)))) * 100;

IRallS = SRa.IRall  * scS;
IRallN = NSRa.IRall * scN;
medS = median(IRallS,4); medN = median(IRallN,4);
qS = prctile(IRallS, [16 84], 4);
qN = prctile(IRallN, [16 84], 4);
writematrix([medS(:,:,mps) qS(:,:,mps,1) qS(:,:,mps,2)], fullfile(out,'ml_adrr_sign.csv'));
writematrix([medN(:,:,mps) qN(:,:,mps,1) qN(:,:,mps,2)], fullfile(out,'ml_adrr_narr.csv'));
writematrix([SRa.accept_rate NSRa.accept_rate scS scN], fullfile(out,'ml_adrr_diag.csv'));

% =====================================================================
% 3. GERTLER & KARADI (2015) — sign + external instrument
% =====================================================================
clearvars -except root out
rng(42,'twister');
f    = fullfile(root,'Replic','GK2015','GK2015_Data.xlsx');
raw  = readcell(f,'Sheet','Sheet1');
mnem = raw(2,2:end);
data = Num2NaN(cellfun(@double, raw(3:end,2:end)));
for ii = 1:numel(mnem); DATA.(mnem{ii}) = data(:,ii); end

ENDO = [DATA.gs1 DATA.logcpi DATA.logip DATA.ebp];
IV   = DATA.ff4_tc;

% Shock 1 is the instrument-identified monetary policy shock; the remaining
% three shocks are matched by sign restrictions (all unrestricted here, so
% the comparison isolates the sign+iv machinery itself).
SIGNiv = zeros(4,3);
SIGNiv(2,1) = -1;   % CPI falls
SIGNiv(3,1) = -1;   % IP falls

VARopt = VARoption;
VARopt.mnem = mnem(1:4);
VARopt.nsteps = 48; VARopt.ndraws = 500; VARopt.sr_hor = 3;
VARopt.pctg = 68;   VARopt.inference = 1; VARopt.sr_draw = 500000;
VARopt.ident = 'sign+iv'; VARopt.R = SIGNiv; VARopt.IV = IV;
SRg = VARmodel(ENDO, 12, 1, VARopt);
fprintf('GK sign+iv done. accept_rate = %.6f\n', SRg.accept_rate);

writematrix([SRg.IRmed(:,:,1) SRg.IRinf(:,:,1) SRg.IRsup(:,:,1)], ...
    fullfile(out,'ml_gk_irf.csv'));
writematrix(SRg.Bmed, fullfile(out,'ml_gk_Bmed.csv'));
writematrix([SRg.accept_rate SRg.ndraws_tried], fullfile(out,'ml_gk_diag.csv'));

% Point IV column for a direct check of fRecoverBIV_cpp
VARopt2 = VARoption; VARopt2.ident = 'iv'; VARopt2.IV = IV;
VARiv = VARmodel(ENDO, 12, 1, VARopt2);
writematrix(VARiv.B(:,1), fullfile(out,'ml_gk_b1.csv'));

fprintf('ALL DONE\n');
