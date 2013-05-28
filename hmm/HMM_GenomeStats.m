 
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

function HMM_GenomeStats

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
global stats

% Initialisations
no_chr = seg_info.no_chr;
no_levels = seg_info.no_levels;
no_genos_ind = seg_info.no_genos_ind;
Tc = seg_info.Tc;
for level_no = 1:no_levels,
    for geno_no = 1:no_genos_ind(level_no),
        a_diag(level_no,geno_no) = params.chi(level_no,geno_no,level_no,geno_no);
    end
end

% Derive LOH levels for each chromosome and genome
for chr_no = 1:no_chr,
    LOH(chr_no,1) = sum(squeeze(sum(alpha{chr_no}(1:no_levels,1,2:Tc(chr_no)+1).*beta{chr_no}(1:no_levels,1,2:Tc(chr_no)+1),1))./c{chr_no}(2:Tc(chr_no)+1)');
end
LOH(no_chr+1,1) = sum(LOH);
LOH = LOH/sum(Tc);

% Derive Complexity values for each chromosome and genome
for chr_no = 1:no_chr,
    comp(chr_no,1) = sum(1 - squeeze(sum(sum(alpha{chr_no}(:,:,1:Tc(chr_no)).*repmat(a_diag,[1,1,Tc(chr_no)]).*b{chr_no}(:,:,1:Tc(chr_no)).*beta{chr_no}(:,:,2:Tc(chr_no)+1),1),2)));
end
comp(no_chr+1,1) = sum(comp);

% Assign output
stats.complexity = comp;
stats.LOH = LOH;

return


