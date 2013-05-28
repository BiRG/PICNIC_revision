
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


function Plot_Genome


% Set Global Variables
global G
global H
global b
global alpha
global beta
global c
global seg_info
global params
global states
global genotypes

% Convert input structure
seg_ends = seg_info.seg_ends;
block_cat = seg_info.block_cat;
block_mean = seg_info.block_mean;
no_blocks = seg_info.no_blocks;
chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;
no_chr = seg_info.no_chr;

% Plot figure
max_height = max(block_mean)*1.5;
figure
hold on;
block_no = 0;
plot([chr_end(no_chr),chr_end(no_chr)],[0,max_height],'k:','LineWidth',0.1);
for chr_no = 1:no_chr,
    plot([chr_start(chr_no),chr_start(chr_no)],[0,max_height],'k:','LineWidth',0.1);
    for seg_no = 1:length(seg_ends{chr_no})-1,
        block_no = block_no + 1;
        if block_cat(block_no) == -1,
            plot([chr_start(chr_no)+seg_ends{chr_no}(seg_no),chr_start(chr_no)+seg_ends{chr_no}(seg_no+1)],[block_mean(block_no),block_mean(block_no)],'y','LineWidth',5);
        elseif  block_cat(block_no) == 1,
            plot([chr_start(chr_no)+seg_ends{chr_no}(seg_no),chr_start(chr_no)+seg_ends{chr_no}(seg_no+1)],[block_mean(block_no),block_mean(block_no)],'m','LineWidth',5);
        elseif  block_cat(block_no) == 2,
            plot([chr_start(chr_no)+seg_ends{chr_no}(seg_no),chr_start(chr_no)+seg_ends{chr_no}(seg_no+1)],[block_mean(block_no),block_mean(block_no)],'r','LineWidth',5);
        elseif  block_cat(block_no) == 3,
            plot([chr_start(chr_no)+seg_ends{chr_no}(seg_no),chr_start(chr_no)+seg_ends{chr_no}(seg_no+1)],[block_mean(block_no),block_mean(block_no)],'k','LineWidth',5);
        else
            plot([chr_start(chr_no)+seg_ends{chr_no}(seg_no),chr_start(chr_no)+seg_ends{chr_no}(seg_no+1)],[block_mean(block_no),block_mean(block_no)],'g','LineWidth',5);
        end    
    end
end
ylim([0 max_height]);

return
