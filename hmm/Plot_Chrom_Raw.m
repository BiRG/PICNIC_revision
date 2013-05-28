function Plot_Chrom_Raw(chr_no)

% Set Global Variables
global G
global H
global seg_info

chr_start = seg_info.chr_start;
chr_end = seg_info.chr_end;

X = seg_info.SNP_pos(chr_start(chr_no):chr_end(chr_no));

figure;
hold on;
plot(X,-2*G{chr_no},'k.','MarkerSize',1);
plot(X,H{chr_no},'k.','MarkerSize',1);
ylim([-2.5 2.5]);


return
