
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


function [init,delta,fn] = Utils_CombFit(peaks2fit,no_points)


m = length(peaks2fit);
n = no_points;
no_steps = 100;
min_peak = min(peaks2fit);
max_peak = max(peaks2fit);
dd_range_mid = (max_peak-min_peak)/(n-1);
min_diff = min([dd_range_mid,peaks2fit(2:m)-peaks2fit(1:m-1)]);
sig = min_diff/15;
rr_range = min_peak-0.5*min_diff:min_diff/no_steps:min_peak+0.5*min_diff;
dd_range = dd_range_mid*0.9:dd_range_mid/no_steps:dd_range_mid*1.1;
init = 0;
delta = 1;
fn = 0;
for rr_val = rr_range,
        
    for dd_val = dd_range,
        fn_val = score_fn(rr_val,dd_val,sig,peaks2fit,m,n);
        if fn_val > fn,
            fn = fn_val;
            init = rr_val;
            delta = dd_val;
        end
    end
end
    

return

function fn_val = score_fn(init,delta,sig,peaks2fit,m,n)


x_mat = repmat(peaks2fit,[n,1]);
ii_mat(1:n,1) = [0:n-1]; 
ii_mat = repmat(ii_mat,[1,m]);
mean_mat = x_mat-(init+ii_mat*delta);

fn_val = prod(sum(sqrt(2*pi)*sig*normpdf(mean_mat,0,sig),1));

return
