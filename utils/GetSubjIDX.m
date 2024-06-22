function [loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ANIMAL_IDs, ANIMAL_VARs] = GetSubjIDX(ROOTDIR, arg_Dir, arg_Meta)
[Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DirectoryAlloc_testedit(ROOTDIR, arg_Dir, arg_Meta);
loopIDX = 1:length(ANIMAL_IDs);
OFCIDX = 1:length(Behavior_files);
subjIDX = cell(1, length(loopIDX));
cnt = 0;
for i = loopIDX(:).'
    cnt = cnt+1;
    subjIDX{cnt} = ANIMAL_VARs.(ANIMAL_IDs{cnt});
end
len = length(OFCIDX);
ANIMAL_IDs
subjIDX
return