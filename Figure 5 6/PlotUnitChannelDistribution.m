[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, 207, 0);
loopIDX = 1:3;
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
OFCIDX = 1:length(Behavior_files);
len = length(OFCIDX);
file_dir = 'C:\Users\NEC_GPU\Downloads\data from Qi\DCZ'; files = natsortfiles(filename_scan(file_dir));
%%
Lumped_ChannelInfo = {};
METAMATRIX = MetaSPK_test_CNO;
cnt = 0;
for i = [subjIDX{3}(end-4:end)]
    cnt = cnt+1;
    load(files{cnt+10})
    goodunit_idx = find(METAMATRIX.Unit_type{i}==1 | METAMATRIX.Unit_type{i}==10 & METAMATRIX.Unit_level{i}>1);
    gathered_ChannelInfo = METAMATRIX.Unit_map{i}(goodunit_idx, 2);
    Lumped_ChannelInfo{end+1, 1} = gathered_ChannelInfo(EncodingCell);
end
LumpedSession_ChannelInfo = cell2mat(Lumped_ChannelInfo);
edges2use = 0:2:200;
% [N, edges] = histcounts(LumpedSession_ChannelInfo, edges2use);
figure; hold on
% plot(edges(1:end-1), N, 'k', 'LineWidth', 1.5)
MEAN_ChannelInfo = [];
for i = 1:length(Lumped_ChannelInfo)
    [MEAN_ChannelInfo(end+1, :), edges] = histcounts(Lumped_ChannelInfo{i}, edges2use, 'Normalization', 'probability');
end
ShadedPlot(edges(1:end-1), nanmean(MEAN_ChannelInfo, 1), 'k', 2, SEM(MEAN_ChannelInfo), [0.8 0.8 0.8])
plot(edges(1:end-1), nanmean(MEAN_ChannelInfo, 1), 'k', 'LineWidth', 1.5)
xlabel('Depth *10 (um)'); ylabel('Fraction')