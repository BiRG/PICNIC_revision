
% Copyright (c) 2008, Genome Research Ltd (GRL).
% All rights reserved.
% Author: Chris Greenman <cdg@sanger.ac.uk>
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above
%       copyright notice, this list of conditions and the following
%       disclaimer in the documentation and/or other materials
%       provided with the distribution.
%     * Neither the name of the <organization> nor the
%       names of its contributors may be used to endorse or promote
%       products derived from this software without specific prior
%       written permission.
%
% THIS SOFTWARE IS PROVIDED BY GRL ``AS IS'' AND ANY EXPRESS OR
% IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED.  IN NO EVENT SHALL GRL BE LIABLE FOR ANY DIRECT,
% INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
% HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
% STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
% OF THE POSSIBILITY OF SUCH DAMAGE.


function Plot_Chrom_HMM(chr_no)


% Set Global Variables
global G
global H
global b
global SNP_ind
global SNP_pos
global alpha
global beta
global c
global seg_info
global params
global states
global genotypes

% Structure conversions
eta = params.eta;
init = params.init;
delta = params.delta;
no_levels = size(params.a0,1);
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
Tc = seg_info.Tc;
cnv = H{chr_no}(1:Tc(chr_no));
geno = 2*no_levels*(G{chr_no}(1:Tc(chr_no))-1)-2;

% Define Geno fitted lines
X_ind = [1:Tc(chr_no)]';
X_SNP = X_ind.*SNP_ind{chr_no};
X_SNP = sort(X_SNP);
X_SNP = X_SNP(end:-1:1);
X_SNP = X_SNP(1:sum(SNP_ind{chr_no}));
X_SNP = X_SNP(end:-1:1);
X_cnv = SNP_pos{chr_no}(1:Tc(chr_no));
X_geno = X_cnv(X_SNP);
Y_geno_top1 = 2/pi*atan(1+(states.cnv{chr_no}-1)*eta);
Y_geno_top2 = 0.5+(1/2-max(0,(states.geno{chr_no}-1)./(states.cnv{chr_no}-1))).*(4/pi*atan(1+(states.cnv{chr_no}-1)*eta)-1);
Y_geno_bot2 = 0.5-(1/2-max(0,(states.geno{chr_no}-1)./(states.cnv{chr_no}-1))).*(4/pi*atan(1+(states.cnv{chr_no}-1)*eta)-1);
Y_geno_bot1 = 1-2/pi*atan(1+(states.cnv{chr_no}-1)*eta);

Y_geno_top1 = 2*no_levels*(Y_geno_top1-1)-2;
Y_geno_top2 = 2*no_levels*(Y_geno_top2-1)-2;
Y_geno_bot2 = 2*no_levels*(Y_geno_bot2-1)-2;
Y_geno_bot1 = 2*no_levels*(Y_geno_bot1-1)-2;

Y_conf = genotypes.geno_state_change{chr_no}*no_levels - no_levels*3 - 2;

% Define Cnv fitted lines
Y_cnv = states.cnv{chr_no}-1;
Y_cnv_low = states.geno{chr_no}-1;

% Get LOH info for distinct colour
LOH_ind = 1*(states.geno{chr_no}==1);
LOH_count = sum(LOH_ind);
XX = sort([1:Tc(chr_no)].*LOH_ind,2,'descend');
X_LOH_ind = XX(1:LOH_count);
X_LOH_ind = X_LOH_ind(end:-1:1);
X_LOH = X_cnv(X_LOH_ind);
Y_cnv_LOH = Y_cnv(X_LOH_ind);
Y_cnv_low_LOH = Y_cnv_low(X_LOH_ind);
Y_geno_top1_LOH = Y_geno_top1(X_LOH_ind);
Y_geno_bot1_LOH = Y_geno_bot1(X_LOH_ind);

% Plot Results
figure;
hold on;
plot(X_cnv,max(-2,(cnv-init)/delta),'k.','MarkerSize',1);
plot(X_geno,geno(X_SNP),'k.','MarkerSize',1);
plot(X_cnv,Y_cnv,'g','MarkerSize',5);
plot(X_cnv,Y_cnv_low,'b','MarkerSize',5);
plot(X_cnv,Y_geno_top1,'k.','MarkerSize',5);
plot(X_cnv,Y_geno_top2,'k.','MarkerSize',5);
plot(X_cnv,Y_geno_bot2,'k.','MarkerSize',5);
plot(X_cnv,Y_geno_bot1,'k.','MarkerSize',5);
plot(X_LOH,Y_geno_top1_LOH,'r.','MarkerSize',5);
plot(X_LOH,Y_geno_bot1_LOH,'r.','MarkerSize',5);
plot(X_cnv,Y_conf,'k');
set(gca,'YTick',0:1.5*no_levels)
title(['SNP info for chr ',num2str(chr_no)])
text(-length(cnv)/30,-2,'A');
text(-length(cnv)/30,-2-2*no_levels,'B');
ylim([-3*no_levels-3 1.5*no_levels]);
xlim([0 SNP_pos{chr_no}(Tc(chr_no))]);

return
