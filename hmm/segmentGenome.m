
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


function segmentGenome


% This function reads in a .feature_intensity found in data/cancer/raw and 1) normalised the data with respect to the normal samples and 2) performs an initial segmentation of the genome

% Input data: This function assumes the directory structure described in the documentation.  The *.feature_intensity files are read in from /data/cancer/raw/*.feature_intensity and places them in /data/cancer/normalised.

%Output.  1) Normalised feature_intensity files go to /data/cancer/normalised
%	  2) Output of segmentation goes to /output/cel_file_name/*

%try
    addpath('../config');
    Features_2CnvGeno;
    HMM_Run_segmentation;
    disp('Finished HMM part A.');
%catch ME
 %   disp('Execution problems.');
 %   return;
%end

return;
