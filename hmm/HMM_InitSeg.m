
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


function HMM_InitSeg


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

% Initialisations
%dbstop if error
no_chr = seg_info.no_chr;
Tc = seg_info.Tc;
levels_max = 10;
lag = 2500;
jump = 50;
t_diff_lim = 0.15;
diff_min = 0.08;
del_cnv = 0.4;

% Find potential breakpoints
for chr_no = 1:no_chr,
    seg_ends{chr_no} = 1;
    inzone_flag = 0;
    for pos = 1+lag:jump:Tc(chr_no)-lag-1,
        X_left = log(H{chr_no}(pos-lag:pos));
        X_right = log(H{chr_no}(pos+1:pos+1+lag));
        mn1 = exp(mean(X_left));
        mn2 = exp(mean(X_right));
        sd1 = std(X_left);
        sd2 = std(X_right);
        t_diff_val = abs((mn1-mn2)/(sd1+sd2));
        diff_val = abs(mn1-mn2);
        if t_diff_val>t_diff_lim & diff_val>diff_min,
            if inzone_flag == 0,
                diff_max = diff_val;
                pos_max = pos;
                inzone_flag = 1;
            elseif diff_val>diff_max & inzone_flag == 1,
                diff_max = diff_val;
                pos_max = pos;
            end
        else
            if inzone_flag > 0,
                inzone_flag = 0;
                seg_ends{chr_no} = [seg_ends{chr_no},pos_max];
            end
        end
    end
    seg_ends{chr_no} = [seg_ends{chr_no},Tc(chr_no)];
end

% Remove breakpoints if segment cnv change is too small
for chr_no = 1:no_chr,
    if length(seg_ends{chr_no})>2,
        seg_shrink_flag = 0;
    else
        seg_shrink_flag = 1;
    end
    while seg_shrink_flag == 0,
        seg_ends2 = 1;
        for seg_no = 2:length(seg_ends{chr_no})-1,
            act_diff = abs(mean(H{chr_no}(seg_ends{chr_no}(seg_no-1):seg_ends{chr_no}(seg_no)))-...
                           mean(H{chr_no}(seg_ends{chr_no}(seg_no):seg_ends{chr_no}(seg_no+1))));
            if act_diff>diff_min,
                seg_ends2 = [seg_ends2,seg_ends{chr_no}(seg_no)];
            end
        end
        seg_ends2 = [seg_ends2,Tc(chr_no)];
        if length(seg_ends2) == length(seg_ends{chr_no}),
            seg_shrink_flag = 1;
        end
        seg_ends{chr_no} = seg_ends2;
    end
end

% Summarise information by segment
block_no = 0;
for chr_no = 1:no_chr,
    XX = [1:Tc(chr_no)]';
    XX = XX.*SNP_ind{chr_no};
    XX = sort(XX);
    XX = XX(end:-1:1);
    XX = XX(1:sum(XX>0));
    XX = XX(end:-1:1);
    for seg_no = 2:length(seg_ends{chr_no}),
        block_no = block_no + 1;
        X1 = sum(XX<seg_ends{chr_no}(seg_no-1))+1;
        X2 = sum(XX<=seg_ends{chr_no}(seg_no));
        ranges = XX(X1:X2);
        [p,q,mu1,mu2,sd1,sd2,fit_score(block_no),block_cat(block_no)] = Utils_ModalQuadFit(G{chr_no}(ranges));
        if chr_no >= no_chr-1,                 % Ignore X,Y chromosomes
            block_cat(block_no) = -1;
        end
        block_mean(block_no) = mean(H{chr_no}(seg_ends{chr_no}(seg_no-1):seg_ends{chr_no}(seg_no)));
        block_std(block_no) = std(H{chr_no}(seg_ends{chr_no}(seg_no-1):seg_ends{chr_no}(seg_no)));
        block_gamma(block_no) = std(log(H{chr_no}(seg_ends{chr_no}(seg_no-1):seg_ends{chr_no}(seg_no))));
        block_size(block_no) = seg_ends{chr_no}(seg_no)-seg_ends{chr_no}(seg_no-1)+1;
        block_mu1(block_no) = mu1;
        block_sd1(block_no) = sd1;
        block_sd2(block_no) = sd2;
        block_q(block_no) = q;
        for bit_no = 1:10,                  % Check block ends for consistent means
            if abs(mean(H{chr_no}(seg_ends{chr_no}(seg_no-1):seg_ends{chr_no}(seg_no-1)+floor(block_size(block_no)/10)))-...
                   mean(H{chr_no}(seg_ends{chr_no}(seg_no)-floor(block_size(block_no)/10):seg_ends{chr_no}(seg_no))))>diff_min,
                block_cat(block_no) = -1;
            end
        end
    end
end
no_blocks = block_no;

% Form structure for output
seg_info.seg_ends = seg_ends;
seg_info.block_cat = block_cat;
seg_info.block_mean = block_mean;
seg_info.block_size = block_size;
seg_info.block_gamma = block_gamma;
seg_info.block_mu1 = block_mu1;
seg_info.block_sd1 = block_sd1;
seg_info.block_sd2 = block_sd2;
seg_info.block_q = block_q;
seg_info.no_blocks = no_blocks;

return

