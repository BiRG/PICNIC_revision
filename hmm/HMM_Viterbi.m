
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


function HMM_Viterbi


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
global p
global params
global states
global genotypes
global segments

% Structure conversion
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
no_chr = seg_info.no_chr;
no_levels = seg_info.no_levels;
no_genos_ind = seg_info.no_genos_ind;
no_states = sum(no_genos_ind);
Tc = seg_info.Tc;
a0 = params.a0;
init = params.init;
delta = params.delta;
gamma = params.gamma;
zeta = params.zeta;
chi = params.chi;
q = params.q;
r = params.r;
eta = params.eta;
sd1 = params.sd1;
sd2 = params.sd2;

% Construct Transition matrix
for level_no1 = 1:no_levels
    for geno_no1 = 1:no_genos_ind(level_no1),
        for level_no2 = 1:no_levels
            for geno_no2 = 1:no_genos_ind(level_no1),
                a(level_no1,geno_no1,level_no2,geno_no2) = chi(level_no1,geno_no1,level_no2,geno_no2);
                a_trans(level_no2,geno_no2,level_no1,geno_no1) = chi(level_no1,geno_no1,level_no2,geno_no2);
            end
        end
    end
end

% Initialisations
eps = 1e-300;

% Loop through chromosomes
for chr_no = 1:no_chr,
    
    chr_no;
    
    % Initialisation
    delta_chr(1:no_levels,1:no_genos_ind(no_levels),1:chr_end(chr_no)-chr_start(chr_no)+1) = -inf;
    phi_level(1:no_levels,1:no_genos_ind(no_levels),1:chr_end(chr_no)-chr_start(chr_no)+1) = 0;
    phi_geno(1:no_levels,1:no_genos_ind(no_levels),1:chr_end(chr_no)-chr_start(chr_no)+1) = 0;
    for level_no = 1:no_levels,
        for geno_no = 1:no_genos_ind(level_no),
            delta_chr(level_no,geno_no,1) = log(eps+a0(level_no,geno_no))+log(eps+b{chr_no}(level_no,geno_no,1));
            phi_level(level_no,geno_no,1) = 0;
            phi_geno(level_no,geno_no,1) = 0;
        end
    end
    
    % Recursion
    for SNP_step = 2:Tc(chr_no),
        for level_no1 = 1:no_levels,
            for geno_no1 = 1:no_genos_ind(level_no1),
                vals_mat = squeeze(delta_chr(:,:,SNP_step-1))+log(eps+squeeze(a(:,:,level_no1,geno_no1)));
                [max_vec,geno_ind] = max(vals_mat,[],2);
                [dum,level_ind_i] = max(max_vec,[],1);
                geno_ind_i = geno_ind(level_ind_i);
                delta_chr(level_no1,geno_no1,SNP_step) = vals_mat(level_ind_i,geno_ind_i) +...
                                                                 log(eps+b{chr_no}(level_no1,geno_no1,SNP_step));
                phi_level(level_no1,geno_no1,SNP_step) = level_ind_i;
                phi_geno(level_no1,geno_no1,SNP_step) = geno_ind_i;
            end
        end
    end
    
    % Termination
    [geno_max,geno_ind] = max(delta_chr(:,:,Tc(chr_no)),[],2);
    [max_val,level_ind] = max(geno_max,[],1);
    geno_ind = geno_ind(level_ind);
    states.cnv{chr_no}(Tc(chr_no)) = level_ind;
    states.geno{chr_no}(Tc(chr_no)) = geno_ind;
    
    % Backtracking
    for SNP_step = Tc(chr_no)-1:-1:1,
        states.cnv{chr_no}(SNP_step) = phi_level(states.cnv{chr_no}(SNP_step+1),states.geno{chr_no}(SNP_step+1),SNP_step+1);
        states.geno{chr_no}(SNP_step) = phi_geno(states.cnv{chr_no}(SNP_step+1),states.geno{chr_no}(SNP_step+1),SNP_step+1);
    end
end

% Output segment info
seg_no = 0;
for chr_no = 1:24,
    pos1 = chr_start(chr_no);
    for t = chr_start(chr_no)+1:chr_end(chr_no),
        t2 = t - chr_start(chr_no) + 1;
        if states.cnv{chr_no}(t2) ~= states.cnv{chr_no}(t2-1) | states.geno{chr_no}(t2) ~= states.geno{chr_no}(t2-1),
            pos2 = t-1;
            seg_no = seg_no + 1;
            segments(seg_no,1:8) = [seg_no,pos1,pos2,chr_no,seg_info.SNP_pos(pos1),seg_info.SNP_pos(pos2),states.geno{chr_no}(t2-1),states.cnv{chr_no}(t2-1)];
            pos1 = t;
        end
        if t == chr_end(chr_no),
            seg_no = seg_no + 1;
            pos2 = t;
            segments(seg_no,1:8) = [seg_no,pos1,pos2,chr_no,seg_info.SNP_pos(pos1),seg_info.SNP_pos(pos2),states.geno{chr_no}(t2),states.cnv{chr_no}(t2)];
        end
    end
end

return
