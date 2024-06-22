%%
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 200, 0);
loopIDX = 1;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end

load([ROOTDIR '2021-22_Attention\NP 2023-12\MetaSPK_test'])
OFCIDX = 1:height(MetaSPK_test);
len = length(OFCIDX);
%%
pairwise_CCF_FR0_Success = cell(len, 2);
pairwise_CCF_FR0_Punish = cell(len, 2);

idx_st = CrossSampling(tt_success, -3);
idx_ed = CrossSampling(tt_success, -1);

for i = 1:len
    [pairwise_CCF_FR0_Success{i, 1}, pairwise_CCF_FR0_Success{i, 2}] = corr(Lumped_PriorSuccess_FR0{i}(:, idx_st:idx_ed)');
    [pairwise_CCF_FR0_Punish{i, 1}, pairwise_CCF_FR0_Punish{i, 2}] = corr(Lumped_PriorPunish_FR0{i}(:, idx_st:idx_ed)');
end

i = 8
figure; hold on
subplot(2, 1, 1);
heatmap(pairwise_CCF_FR0_Success{i, 1}); colormap(jet); caxis([-0.8 1]); grid off
subplot(2, 1, 2);
heatmap(pairwise_CCF_FR0_Punish{i, 1}); colormap(jet); caxis([-0.8 1]); grid off
%%
% figure; hold on
% for i = 1:len
%     subplot(4, 4, i); hold on
%     pvalue_mask = pairwise_CCF_FR0_Success{i, 2};
%     pvalue_mask(find(pvalue_mask>0.05)) = NaN;
%     flattened_CCF = pairwise_CCF_FR0_Success{i, 1}.*pvalue_mask;
%     flattened_CCF = flattened_CCF(:);
%     [N, edges] = histcounts(flattened_CCF, 'Normalization','probability');
%     edges = edges(2:end) - (edges(2)-edges(1))/2;
%     plot(edges, N, 'g');
%     
%     pvalue_mask = pairwise_CCF_FR0_Punish{i, 2};
%     pvalue_mask(find(pvalue_mask>0.05)) = NaN;
%     flattened_CCF = pairwise_CCF_FR0_Punish{i, 1}.*pvalue_mask;
%     flattened_CCF = flattened_CCF(:);
%     [N, edges] = histcounts(flattened_CCF, 'Normalization','probability');
%     edges = edges(2:end) - (edges(2)-edges(1))/2;
%     plot(edges, N, 'm');
% end

figure; hold on
for i = 1:len
%     subplot(4, 4, i); hold on
    pvalue_mask = triu(pairwise_CCF_FR0_Success{i, 2}, 1);
    pvalue_mask(find(pvalue_mask==0 | pvalue_mask>0.05)) = NaN;
    pvalue_mask(find(pvalue_mask<=0.05)) = 1;
    significant_CCF_Success = pairwise_CCF_FR0_Success{i, 1}.*pvalue_mask;
    
    pvalue_mask = triu(pairwise_CCF_FR0_Punish{i, 2}, 1);
    pvalue_mask(find(pvalue_mask==0 | pvalue_mask>0.05)) = NaN;
    pvalue_mask(find(pvalue_mask<=0.05)) = 1;
    significant_CCF_Punish = pairwise_CCF_FR0_Punish{i, 1}.*pvalue_mask;
    scatter(nanmean(abs(significant_CCF_Punish(:))), nanmean(abs(significant_CCF_Success(:))), 20, 'k', 'o')
end
% for i = 1:len
%     subplot(4, 4, i); hold on
%     scatter(
    
    
    
    

