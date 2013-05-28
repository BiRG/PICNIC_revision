
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


function HMM_InitParamSet


% Set Global Variables
global G
global H
global b
global SNP_ind
global alpha
global beta
global c
global seg_info
global params
global states
global genotypes
global del_in
global del_out

% Convert input structure
seg_ends = seg_info.seg_ends;
block_cat = seg_info.block_cat;
block_mean = seg_info.block_mean;
block_size = seg_info.block_size;
block_mu1 = seg_info.block_mu1;
no_blocks = seg_info.no_blocks;
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
levels = seg_info.levels;
no_levels = seg_info.no_levels;
no_genos_ind = seg_info.no_genos_ind;
no_genos_max = seg_info.no_genos_max;
no_SNPs = seg_info.no_SNPs;
chr_len = seg_info.Tc;

% Initialisations
init = levels(1);
delta = levels(2)-levels(1);

% Get Initial distribution 'vector' (really a matrix)
pii(1:no_levels) = 1;
qii(1:no_levels,1:no_genos_max) = 0;
for block_no = 1:no_blocks,
    [dumm,ind2use] = min(abs(block_mean(block_no)-levels));
    block_level(block_no) = levels(ind2use);
    block_ind(block_no) = ind2use;
    pii(ind2use) = pii(ind2use) + block_size(block_no);
    if block_cat(block_no) == 1,
        qii(ind2use,1) = qii(ind2use,1) + block_size(block_no);
    end
end
qii = qii./repmat(pii'+eps,[1,no_genos_max]);
qii(1:2,1) = 1;
for level_no = 1:no_levels,
    for geno_no = 1:no_genos_ind(level_no),
        qii(level_no,geno_no) = pii(level_no)/no_genos_ind(level_no);
    end
end
state_count = pii;
pii = pii'/sum(pii);
init_matrix = qii/sum(sum(qii));

% Get transition 'matrix' parameter (really a tensor)
chi_change = 1e-18;
chi_change2del = 1e-14;
chi_change4del = 1e-8;
for level_no1 = 1:no_levels,
    for geno_no1 = 1:no_genos_ind(level_no1),
        for level_no2 = 1:no_levels,
            for geno_no2 = 1:no_genos_ind(level_no2),
                if level_no1==level_no2 & geno_no1==geno_no2,
                    chi(level_no1,geno_no1,level_no2,geno_no2) = 1;
                else
                    if level_no1==1,
                        chi(level_no1,geno_no1,level_no2,geno_no2) = chi_change4del;
                    else
                        if level_no2==1,
                            chi(level_no1,geno_no1,level_no2,geno_no2) = chi_change2del;
                        else
                            chi(level_no1,geno_no1,level_no2,geno_no2) = chi_change;
                        end
                    end
                end
            end
        end
        chi(level_no1,geno_no1,:,:) = chi(level_no1,geno_no1,:,:)/sum(sum(chi(level_no1,geno_no1,:,:)));
    end
end

% Emission Parameters
gamma = 0.3;
sd1 = 0.01;
sd2 = 0.04;
q = 0.3;
eta = 0.7;
for block_no = 1:no_blocks,
    level2use = 1+(block_mean(block_no)-init)/delta;
    if level2use < 6 & level2use > 1,
        eta = [eta,(tan(pi/2*(block_mu1(block_no)+1/2))-1)/max(1,(level2use-1))];
        q = [q,seg_info.block_q(block_no)];
        sd1 = [sd1,seg_info.block_sd1(block_no)];
        sd2 = [sd2,seg_info.block_sd2(block_no)];
        gamma = [gamma,seg_info.block_gamma(block_no)];
    end
end

%Set output structure
params.a0 = init_matrix;
params.chi = chi;
params.init = init;
params.delta = delta;
params.gamma = median(gamma);
params.zeta = 10;
params.q = median(q);
params.r = 0.995;
params.eta = median(eta);
params.sd1 = median(sd1);
params.sd2 = median(sd2);

return

