[e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 7, 0, -1, 2);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 14, 0, 0, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 6, 0, 1, 1);
% [e1, loopIDX, OFCIDX, subjIDX, len, Behavior_files, Phot_files, Pupil_files, ...
%     MetaData_files, MetaDataX_files, ANIMAL_IDs, ANIMAL_VARs] = DoMeFavor(ROOTDIR, 13, 0, 1, 1);
session_Pupil = cell(len, 1);
session_Pupil_RC1 = cell(len, 1);
session_Pupil_RC2 = cell(len, 1);
session_Pupil_RC3 = cell(len, 1);

session_bdPupil = cell(len, 1);

session_phiRC3 = cell(len, 1);

session_dPupil = cell(len, 1);
session_dPupil_RC = cell(len, 1);

session_HiPupil = cell(len, 1);
session_HiPupilRC1 = cell(len, 1);
session_HiPupilRC2 = cell(len, 1);
session_HiPupilRC3 = cell(len, 1);
session_phiHiPupil = cell(len, 1);
session_HiPupil_LDA_acc = NaN(len, 1);


session_LoPupil = cell(len, 1);
session_LoPupilRC1 = cell(len, 1);
session_LoPupilRC2 = cell(len, 1);
session_LoPupilRC3 = cell(len, 1);

session_bdPupil_atOutcome_P = cell(len, 1);
session_bdPupil_atOutcome_R = cell(len, 1);

session_phiPupil = cell(len, 1);

session_dilaPupil = cell(len, 1);
session_consPupil = cell(len, 1);

% session_Pupil_atTone_P = cell(len, 1); 
% session_Pupil_atTone_R = cell(len, 1); 
% session_dPupil_atTone_P = cell(len, 1);
% session_dPupil_atTone_R = cell(len, 1);
% session_HiPupil_atTone_P = cell(len, 1); 
% session_HiPupil_atTone_R = cell(len, 1);
% session_LoPupil_atTone_P = cell(len, 1); 
% session_LoPupil_atTone_R = cell(len, 1);
session_Pupil_atOutcome_NS = cell(len, 1);
session_dPupil_atOutcome_NS = cell(len, 1);
session_Pupil_Free_NS = cell(len, 1);
session_dPupil_Free_NS = cell(len, 1);

session_Pupil_atOutcome_P = cell(len, 1); 
session_PupilRC1_atOutcome_P = cell(len, 1); 
session_PupilRC2_atOutcome_P = cell(len, 1);
session_PupilRC3_atOutcome_P = cell(len, 1);
session_Pupil_atOutcome_R = cell(len, 1); 
session_PupilRC1_atOutcome_R = cell(len, 1); 
session_PupilRC2_atOutcome_R = cell(len, 1); 
session_PupilRC3_atOutcome_R = cell(len, 1); 

session_dPupil_atOutcome_P = cell(len, 1);
session_dPupil_atOutcome_R = cell(len, 1);

session_HiPupil_atOutcome_P = cell(len, 1); 
session_HiPupilRC1_atOutcome_P = cell(len, 1);
session_HiPupilRC2_atOutcome_P = cell(len, 1);
session_HiPupilRC3_atOutcome_P = cell(len, 1);
session_HiPupil_atOutcome_R = cell(len, 1);
session_HiPupilRC1_atOutcome_R = cell(len, 1);
session_HiPupilRC2_atOutcome_R = cell(len, 1);
session_HiPupilRC3_atOutcome_R = cell(len, 1);

session_LoPupil_atOutcome_P = cell(len, 1); 
session_LoPupilRC1_atOutcome_P = cell(len, 1);
session_LoPupilRC2_atOutcome_P = cell(len, 1);
session_LoPupilRC3_atOutcome_P = cell(len, 1);
session_LoPupil_atOutcome_R = cell(len, 1);
session_LoPupilRC1_atOutcome_R = cell(len, 1);
session_LoPupilRC2_atOutcome_R = cell(len, 1);
session_LoPupilRC3_atOutcome_R = cell(len, 1);

session_phiPupil_atOutcome_P = cell(len, 1);
session_phiPupil_atOutcome_R = cell(len, 1);

session_dilaPupil_P = cell(len, 1);
session_consPupil_P = cell(len, 1);
session_dilaPupil_R = cell(len, 1);
session_consPupil_R = cell(len, 1);

session_pupilIDX_Last5sOC_P = cell(len, 1);
session_pupilIDX_Last5sOC_R = cell(len, 1);
session_pupilIDX_TonePeriod_P = cell(len, 1);
session_pupilIDX_TonePeriod_R = cell(len, 1);
session_pupilIDX_LickFreePeriod_P = cell(len, 1);
session_pupilIDX_LickFreePeriod_R = cell(len, 1);
session_pupilIDX_Trial2Outcome_P = cell(len, 1);
session_pupilIDX_Trial2Outcome_R = cell(len, 1);
session_pupilIDX_shuffle = cell(len, 1);

bandPupil = 0.1;
band = [0.1 1];
winsz = 5;
tic
OPPTWIN = {'Outcome05', 'Outcome075', 'Outcome1', 'Outcome15', 'Outcome2'};
option = 3; 
bandSensor = [0.4 0.8];
WINspkcnt_st = -3;
WINspkcnt_ed = -1;
WINfree_st = -5; 
WINfree_ed = 0;
EDGE_FR = -3.5:0.5:-1;

cnt = 0;

for i = OFCIDX(:).'
    cnt = cnt+1; 
    r_pupil = e1.MetaData.Pupil_size{i};
    z_pupil = zscore(r_pupil);
    z_dpupil = zscore(diff(r_pupil));
    
    
    z_pupil_notnan = zscore(Filter(r_pupil(~isnan(r_pupil)), 10, 2, bandPupil, 'high'));
    z_Hipupil = r_pupil; z_Hipupil(~isnan(r_pupil)) = z_pupil_notnan;
    z_pupil_notnan = zscore(Filter(r_pupil(~isnan(r_pupil)), 10, 2, bandPupil, 'low'));
    z_Lopupil = r_pupil; z_Lopupil(~isnan(r_pupil)) = z_pupil_notnan;
    
    z_pupil_notnan = zscore(Filter(r_pupil(~isnan(r_pupil)), 10, 2, band, 'bandpass'));
    z_bdpupil = r_pupil; z_bdpupil(~isnan(r_pupil)) = z_pupil_notnan;
    phipupil_notnan = angle(hilbert(z_pupil_notnan));
    phipupil = r_pupil; phipupil(~isnan(r_pupil)) = phipupil_notnan;

    timestamp = e1.MetaData.Pupil_time{i};
    mx = load(MetaDataX_files{i});
    if isempty(r_pupil) || ~isempty(find(isnan(r_pupil))) || ismember(i, [2 3 8 10 129])
        session_Pupil{cnt} = [];
        session_Pupil_RC1{cnt} = [];
        session_Pupil_RC2{cnt} = [];
        session_Pupil_RC3{cnt} = [];
        session_bdPupil{cnt} = [];
        session_phiRC3{cnt} = [];
        session_dPupil{cnt} = [];
        session_HiPupil{cnt} = [];
        session_HiPupilRC1{cnt} = [];
        session_HiPupilRC2{cnt} = [];
        session_HiPupilRC3{cnt} = [];
        session_phiHiPupil{cnt} = [];
%         session_HiPupil_LDA_acc(cnt) = [];
        session_LoPupil{cnt} = [];
        session_LoPupilRC1{cnt} = cell(len, 1);
        session_LoPupilRC2{cnt} = cell(len, 1);
        session_LoPupilRC3{cnt} = cell(len, 1);
        session_phiPupil{cnt} = [];
        session_dilaPupil{cnt} = [];
        session_consPupil{cnt} = [];
        session_Pupil_atOutcome_NS{cnt} = [];
        session_dPupil_atOutcome_NS{cnt} = [];
        session_Pupil_Free_NS{cnt} = [];
        session_dPupil_Free_NS{cnt} = [];
        session_Pupil_atOutcome_P{cnt} = [];
        session_Pupil_atOutcome_R{cnt} = [];
        session_dPupil_atOutcome_P{cnt} = [];
        session_dPupil_atOutcome_R{cnt} = [];
        session_HiPupil_atOutcome_P{cnt} = [];
        session_HiPupilRC1_atOutcome_P{cnt} = [];
        session_HiPupilRC2_atOutcome_P{cnt} = [];
        session_HiPupilRC3_atOutcome_P{cnt} = [];
        session_HiPupil_atOutcome_R{cnt} = [];
        session_HiPupilRC1_atOutcome_R{cnt} = [];
        session_HiPupilRC2_atOutcome_R{cnt} = [];
        session_HiPupilRC3_atOutcome_R{cnt} = [];
        session_LoPupil_atOutcome_P{cnt} = [];
        session_LoPupil_atOutcome_R{cnt} = [];
        session_pupilIDX_Last5sOC_P{cnt} = [];
        session_pupilIDX_Last5sOC_R{cnt} = [];
        session_pupilIDX_TonePeriod_P{cnt} = {};
        session_pupilIDX_TonePeriod_R{cnt} = {};
        session_pupilIDX_LickFreePeriod_P{cnt} = [];
        session_pupilIDX_LickFreePeriod_R{cnt} = [];
        session_pupilIDX_Trial2Outcome_P{cnt} = {};
        session_pupilIDX_Trial2Outcome_R{cnt} = {};
        session_pupilIDX_shuffle{cnt} = [];
    else
        [RC, PC, LAMBDA, RHO] = SSA(z_pupil, 15);
        RC1 = RC(:, 1); RC2 = RC(:, 2); RC3 = RC(:, 3);
        session_Pupil{cnt} = z_pupil;
        session_Pupil_RC1{cnt} = RC1;
        session_Pupil_RC2{cnt} = RC2;
        session_Pupil_RC3{cnt} = RC3;
        session_bdPupil{cnt} = z_bdpupil;
%         session_phiRC3{cnt} = angle(hilbert(RC3));
        session_dPupil{cnt} = z_dpupil;
%         session_dPupil_RC{cnt} = diff(RC1);
%         [RC, PC, LAMBDA, RHO] = SSA(z_Hipupil, 15);
%         RC1 = RC(:, 1); RC2 = RC(:, 2); RC3 = RC(:, 3);
        session_HiPupil{cnt} = z_Hipupil;
%         session_HiPupilRC1{cnt} = RC1;
%         session_HiPupilRC2{cnt} = RC2;
%         session_HiPupilRC3{cnt} = RC3;
        session_phiHiPupil{cnt} = phipupil;
        
%         [RC, PC, LAMBDA, RHO] = SSA(z_Lopupil, 15);
%         RC1 = RC(:, 1); RC2 = RC(:, 2); RC3 = RC(:, 3);
        session_LoPupil{cnt} = z_Lopupil;
%         session_LoPupilRC1{cnt} = RC1;
%         session_LoPupilRC2{cnt} = RC2;
%         session_LoPupilRC3{cnt} = RC3;
        session_phiPupil{cnt} = phipupil;
        [dila, cons] = CalcCycles(z_bdpupil);
        session_dilaPupil{cnt} = dila;
        session_consPupil{cnt} = cons;
        
%         for sf = 1:200 
%             rng('shuffle')                    % shuffle for pupil signal part
%             session_pupilIDX_shuffle{cnt}(end+1, :) = randperm(length(timestamp));
%         end
        for k = 2:height(mx.MetaDataX)
            stTone_pupil = CrossSampling(timestamp, mx.MetaDataX.Tone_Onset(k));
            stTrlOnset_pupil = CrossSampling(timestamp, mx.MetaDataX.Trial_Onset(k));
            if mx.MetaDataX.Punish_Onset(k)-mx.MetaDataX.Tone_Onset(k)>5
                stOutcome_pupil = CrossSampling(timestamp, mx.MetaDataX.Punish_Onset(k));
                session_Pupil_atOutcome_P{cnt}(end+1, :) = z_pupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_PupilRC1_atOutcome_P{cnt}(end+1, :) = RC1(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_PupilRC2_atOutcome_P{cnt}(end+1, :) = RC2(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_PupilRC3_atOutcome_P{cnt}(end+1, :) = RC3(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_dPupil_atOutcome_P{cnt}(end+1, :) = z_dpupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_HiPupil_atOutcome_P{cnt}(end+1, :) = z_Hipupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
%                 session_HiPupilRC1_atOutcome_P{cnt}(end+1, :) = RC1(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
%                 session_HiPupilRC2_atOutcome_P{cnt}(end+1, :) = RC2(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
%                 session_HiPupilRC3_atOutcome_P{cnt}(end+1, :) = RC3(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_LoPupil_atOutcome_P{cnt}(end+1, :) = z_Lopupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_bdPupil_atOutcome_P{cnt}(end+1, :) = z_bdpupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_phiPupil_atOutcome_P{cnt}(end+1, :) = phipupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_Pupil_atOutcome_NS{cnt}(end+1, :) = z_pupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_dPupil_atOutcome_NS{cnt}(end+1, :) = z_dpupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_Pupil_Free_NS{cnt}(end+1, :) = z_pupil(stTone_pupil-10*10:stTone_pupil+1*10)';
                session_dPupil_Free_NS{cnt}(end+1, :) = z_dpupil(stTone_pupil-10*10:stTone_pupil+1*10)';
                idx = intersect(dila, stOutcome_pupil-5*10:stOutcome_pupil);
                session_dilaPupil_P{cnt} = [session_dilaPupil_P{cnt}; idx];
                idx = intersect(cons, stOutcome_pupil-5*10:stOutcome_pupil);
                session_consPupil_P{cnt} = [session_consPupil_P{cnt}; idx];
                
                session_pupilIDX_Last5sOC_P{cnt}(end+1, :) = stOutcome_pupil-50:stOutcome_pupil+abs(tmin)/100+abs(tmax)/100;
                session_pupilIDX_TonePeriod_P{cnt}{end+1, 1} = stTone_pupil:stOutcome_pupil+abs(tmin)/100+abs(tmax)/100;
                session_pupilIDX_LickFreePeriod_P{cnt}{end+1, 1} = stTrlOnset_pupil+WINfree_st*10:stTone_pupil+WINfree_ed*10+abs(tmin)/100+abs(tmax)/100;
                session_pupilIDX_Trial2Outcome_P{cnt}{end+1, 1} = stTrlOnset_pupil:stOutcome_pupil+abs(tmin)/100+abs(tmax)/100;
            end
            if mx.MetaDataX.(OPPTWIN{option})(k)==1
                stOutcome_pupil = CrossSampling(timestamp, mx.MetaDataX.Reward_Onset(k));
                session_Pupil_atOutcome_R{cnt}(end+1, :) = z_pupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_PupilRC1_atOutcome_R{cnt}(end+1, :) = RC1(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_PupilRC2_atOutcome_R{cnt}(end+1, :) = RC2(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_PupilRC3_atOutcome_R{cnt}(end+1, :) = RC3(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_dPupil_atOutcome_R{cnt}(end+1, :) = z_dpupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_HiPupil_atOutcome_R{cnt}(end+1, :) = z_Hipupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
%                 session_HiPupilRC1_atOutcome_R{cnt}(end+1, :) = RC1(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
%                 session_HiPupilRC2_atOutcome_R{cnt}(end+1, :) = RC2(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
%                 session_HiPupilRC3_atOutcome_R{cnt}(end+1, :) = RC3(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_LoPupil_atOutcome_R{cnt}(end+1, :) = z_Lopupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_bdPupil_atOutcome_R{cnt}(end+1, :) = z_bdpupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_phiPupil_atOutcome_R{cnt}(end+1, :) = phipupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_Pupil_atOutcome_NS{cnt}(end+1, :) = z_pupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_dPupil_atOutcome_NS{cnt}(end+1, :) = z_dpupil(stOutcome_pupil-5*10:stOutcome_pupil+1*10)';
                session_Pupil_Free_NS{cnt}(end+1, :) = z_pupil(stTone_pupil-10*10:stTone_pupil+1*10)';
                session_dPupil_Free_NS{cnt}(end+1, :) = z_dpupil(stTone_pupil-10*10:stTone_pupil+1*10)';
                idx = intersect(dila, stOutcome_pupil-5*10:stOutcome_pupil);
                session_dilaPupil_R{cnt} = [session_dilaPupil_R{cnt}; idx];
                idx = intersect(cons, stOutcome_pupil-5*10:stOutcome_pupil);
                session_consPupil_R{cnt} = [session_consPupil_R{cnt}; idx];
                
                session_pupilIDX_Last5sOC_R{cnt}(end+1, :) = stOutcome_pupil-50:stOutcome_pupil+abs(tmin)/100+abs(tmax)/100;
                session_pupilIDX_TonePeriod_R{cnt}{end+1, 1} = stTone_pupil:stOutcome_pupil+abs(tmin)/100+abs(tmax)/100;
                session_pupilIDX_LickFreePeriod_R{cnt}{end+1, 1} = stTrlOnset_pupil+WINfree_st*10:stTone_pupil+WINfree_ed*10+abs(tmin)/100+abs(tmax)/100;
                session_pupilIDX_Trial2Outcome_R{cnt}{end+1, 1} = stTone_pupil:stOutcome_pupil+abs(tmin)/100+abs(tmax)/100;
            end
        end
%         num0 = size(session_HiPupil_atOutcome_P{cnt}, 1);
%         num1 = size(session_HiPupil_atOutcome_R{cnt}, 1);
%         raw_data = [session_HiPupil_atOutcome_P{cnt}(:, 1:40); session_HiPupil_atOutcome_R{cnt}(:, 1:40)];
%         test_labels = [zeros(num0, 1); ones(num1, 1)];
%         session_HiPupil_LDA_acc(cnt) = LDA(raw_data, test_labels);
    end
end
toc
%%
Lumped_Pupil_atOutcome_P = NaN(length(OFCIDX), 61); % Punished but withheld >4s NE
Lumped_PupilRC1_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_PupilRC2_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_PupilRC3_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_Pupil_atOutcome_R = NaN(length(OFCIDX), 61);
Lumped_PupilRC1_atOutcome_R = NaN(length(OFCIDX), 61);
Lumped_PupilRC2_atOutcome_R = NaN(length(OFCIDX), 61);
Lumped_PupilRC3_atOutcome_R = NaN(length(OFCIDX), 61);

Lumped_dPupil_atOutcome_P = NaN(length(OFCIDX), 61); % Punished but withheld >4s Ach
Lumped_dPupil_atOutcome_R = NaN(length(OFCIDX), 61);

Lumped_HiPupil_atOutcome_P = NaN(length(OFCIDX), 61); % Rewarded NE
Lumped_HiPupilRC1_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_HiPupilRC2_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_HiPupilRC3_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_HiPupil_atOutcome_R = NaN(length(OFCIDX), 61);
Lumped_HiPupilRC1_atOutcome_R = NaN(length(OFCIDX), 61);
Lumped_HiPupilRC2_atOutcome_R = NaN(length(OFCIDX), 61);
Lumped_HiPupilRC3_atOutcome_R = NaN(length(OFCIDX), 61);

Lumped_LoPupil_atOutcome_P = NaN(length(OFCIDX), 61); % Rewarded Ach
Lumped_LoPupil_atOutcome_R = NaN(length(OFCIDX), 61);

Lumped_bdPupil_atOutcome_P = NaN(length(OFCIDX), 61);
Lumped_bdPupil_atOutcome_R = NaN(length(OFCIDX), 61);

for i = 1:length(OFCIDX)
    if isempty(session_Pupil_atOutcome_P{i}) || isempty(session_Pupil_atOutcome_R{i}) || isempty(session_dPupil_atOutcome_P{i}) || isempty(session_dPupil_atOutcome_R{i})...
             || isempty(session_HiPupil_atOutcome_P{i}) || isempty(session_HiPupil_atOutcome_R{i}) || isempty(session_LoPupil_atOutcome_P{i}) || isempty(session_LoPupil_atOutcome_R{i})
        continue;
    end
    Lumped_Pupil_atOutcome_P(i, :) = nanmean(session_Pupil_atOutcome_P{i}, 1);
    Lumped_PupilRC1_atOutcome_P(i, :) = nanmean(session_PupilRC1_atOutcome_P{i}, 1);
    Lumped_PupilRC2_atOutcome_P(i, :) = nanmean(session_PupilRC2_atOutcome_P{i}, 1);
    Lumped_PupilRC3_atOutcome_P(i, :) = nanmean(session_PupilRC3_atOutcome_P{i}, 1);
    Lumped_Pupil_atOutcome_R(i, :) = nanmean(session_Pupil_atOutcome_R{i}, 1);
    Lumped_PupilRC1_atOutcome_R(i, :) = nanmean(session_PupilRC1_atOutcome_R{i}, 1);
    Lumped_PupilRC2_atOutcome_R(i, :) = nanmean(session_PupilRC2_atOutcome_R{i}, 1);
    Lumped_PupilRC3_atOutcome_R(i, :) = nanmean(session_PupilRC3_atOutcome_R{i}, 1);
    Lumped_dPupil_atOutcome_P(i, :) = nanmean(session_dPupil_atOutcome_P{i}, 1);
    Lumped_dPupil_atOutcome_R(i, :) = nanmean(session_dPupil_atOutcome_R{i}, 1);
    Lumped_HiPupil_atOutcome_P(i, :) = nanmean(session_HiPupil_atOutcome_P{i}, 1);
%     Lumped_HiPupilRC1_atOutcome_P(i, :) = nanmean(session_HiPupilRC1_atOutcome_P{i}, 1);
%     Lumped_HiPupilRC2_atOutcome_P(i, :) = nanmean(session_HiPupilRC2_atOutcome_P{i}, 1);
%     Lumped_HiPupilRC3_atOutcome_P(i, :) = nanmean(session_HiPupilRC3_atOutcome_P{i}, 1);
    Lumped_HiPupil_atOutcome_R(i, :) = nanmean(session_HiPupil_atOutcome_R{i}, 1);
%     Lumped_HiPupilRC1_atOutcome_R(i, :) = nanmean(session_HiPupilRC1_atOutcome_R{i}, 1);
%     Lumped_HiPupilRC2_atOutcome_R(i, :) = nanmean(session_HiPupilRC2_atOutcome_R{i}, 1);
%     Lumped_HiPupilRC3_atOutcome_R(i, :) = nanmean(session_HiPupilRC3_atOutcome_R{i}, 1);
    Lumped_LoPupil_atOutcome_P(i, :) = nanmean(session_LoPupil_atOutcome_P{i}, 1);
    Lumped_LoPupil_atOutcome_R(i, :) = nanmean(session_LoPupil_atOutcome_R{i}, 1);
    Lumped_bdPupil_atOutcome_P(i, :) = nanmean(session_bdPupil_atOutcome_P{i}, 1);
    Lumped_bdPupil_atOutcome_R(i, :) = nanmean(session_bdPupil_atOutcome_R{i}, 1);
end
%%
TT2 = (-5*10:1*10)/10;
figure; hold on
% subplot(4, 1, 1); hold on
ShadedPlot(TT2, nanmean(Lumped_HiPupil_atOutcome_P, 1), [0 0 1], 1, SEM(Lumped_HiPupil_atOutcome_P), [0.8 0.8 0.8])
plot(TT2, nanmean(Lumped_HiPupil_atOutcome_P, 1), '--k', 'LineWidth', 1)
ShadedPlot(TT2, nanmean(Lumped_HiPupil_atOutcome_R, 1), [0 0 1], 1, SEM(Lumped_HiPupil_atOutcome_R), [0.8 0.8 0.8])
plot(TT2, nanmean(Lumped_HiPupil_atOutcome_R, 1), '-k', 'LineWidth', 1)
ylabel('Z-score'); vline(0, '-k'); title('HPF pupil')
% subplot(4, 1, 2); hold on
% ShadedPlot(TT2, nanmean(Lumped_HiPupilRC1_atOutcome_P, 1), [1 0 0], 1, SEM(Lumped_HiPupilRC1_atOutcome_P), [0.8 0.8 0.8])
% plot(TT2, nanmean(Lumped_HiPupilRC1_atOutcome_P, 1), '--k', 'LineWidth', 1)
% ShadedPlot(TT2, nanmean(Lumped_HiPupilRC1_atOutcome_R, 1), [1 0 0], 1, SEM(Lumped_HiPupilRC1_atOutcome_R), [0.8 0.8 0.8])
% plot(TT2, nanmean(Lumped_HiPupilRC1_atOutcome_R, 1), '-k', 'LineWidth', 1)
% ylabel('Z-score'); vline(0, '-k'); title('HPF pupil PC1')
% subplot(4, 1, 3); hold on
% ShadedPlot(TT2, nanmean(Lumped_HiPupilRC2_atOutcome_P, 1), [1 0 0], 1, SEM(Lumped_HiPupilRC2_atOutcome_P), [0.8 0.8 0.8])
% plot(TT2, nanmean(Lumped_HiPupilRC2_atOutcome_P, 1), '--k', 'LineWidth', 1)
% ShadedPlot(TT2, nanmean(Lumped_HiPupilRC2_atOutcome_R, 1), [1 0 0], 1, SEM(Lumped_HiPupilRC2_atOutcome_R), [0.8 0.8 0.8])
% plot(TT2, nanmean(Lumped_HiPupilRC2_atOutcome_R, 1), '-k', 'LineWidth', 1)
% ylabel('Z-score'); vline(0, '-k'); title('HPF pupil PC2')
% subplot(4, 1, 4); hold on
% ShadedPlot(TT2, nanmean(Lumped_HiPupilRC3_atOutcome_P, 1), [1 0 0], 1, SEM(Lumped_HiPupilRC3_atOutcome_P), [0.8 0.8 0.8])
% plot(TT2, nanmean(Lumped_HiPupilRC3_atOutcome_P, 1), '--k', 'LineWidth', 1)
% ShadedPlot(TT2, nanmean(Lumped_HiPupilRC3_atOutcome_R, 1), [1 0 0], 1, SEM(Lumped_HiPupilRC3_atOutcome_R), [0.8 0.8 0.8])
% plot(TT2, nanmean(Lumped_HiPupilRC3_atOutcome_R, 1), '-k', 'LineWidth', 1)
% xlabel('Time after Outcome (s)'); ylabel('Z-score'); vline(0, '-k'); title('HPF pupil PC3')


%% plot a specific session's pupil size traces
figure; hold on
s = 139;
plot((1:size(Lumped_HiPupil_atOutcome_P, 2))/10-5, Lumped_HiPupil_atOutcome_P(s, :), '--k', 'LineWidth', 1)
plot((1:size(Lumped_HiPupil_atOutcome_P, 2))/10-5, Lumped_HiPupil_atOutcome_R(s, :), '-k', 'LineWidth', 1)

%%
figure; hold on
WINst = 2*10+1;
WINed = 4*10;
MEAN_pupil_P = [];
MEAN_pupil_R = [];
for subj = 1:length(subjIDX)
    MEAN_pupil_P(end+1, 1) = nanmean(nanmean(Lumped_HiPupil_atOutcome_P(subjIDX{subj}, WINst:WINed), 2));
    MEAN_pupil_R(end+1, 1) = nanmean(nanmean(Lumped_HiPupil_atOutcome_R(subjIDX{subj}, WINst:WINed), 2));
end

bar(1, nanmean(MEAN_pupil_P))
bar(2, nanmean(MEAN_pupil_R))
for i = 1:length(subjIDX)
    plot(1:2, [MEAN_pupil_P(i) MEAN_pupil_R(i)], '-ok')
end
errorbar(1:2, [nanmean(MEAN_pupil_P) nanmean(MEAN_pupil_R)], [SEM(MEAN_pupil_P) SEM(MEAN_pupil_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Mean (z-score)')

%%
figure; hold on
WINst1 = 2.5*10;
WINed1 = 3.5*10;
WINst2 = 4.5*10;
WINed2 = 5.5*10;
Slope_P = [];
Slope_R = [];
for i = 1:length(subjIDX)
    trace_P = medfilt1(nanmean(Lumped_HiPupil_atOutcome_P(subjIDX{i}, :), 1), 3);
    Slope_P(end+1, 1) = (nanmean(trace_P(WINst1:WINed1))-nanmean(trace_P(WINst2:WINed2)))/((WINed1-WINed2)/10);
    trace_R = medfilt1(nanmean(Lumped_HiPupil_atOutcome_R(subjIDX{i}, :), 1), 3);
    Slope_R(end+1, 1) = (nanmean(trace_R(WINst1:WINed1))-nanmean(trace_R(WINst2:WINed2)))/((WINed1-WINed2)/10);
end

bar(1, nanmean(Slope_P))
bar(2, nanmean(Slope_R))
for i = 1:length(subjIDX)
    plot(1:2, [Slope_P(i) Slope_R(i)], '-ok')
end
errorbar(1:2, [nanmean(Slope_P) nanmean(Slope_R)], [SEM(Slope_P) SEM(Slope_R)], 'm', 'LineStyle', 'none', 'CapSize', 6, 'LineWidth', 1);
ylabel('Slope (a.u.)')

























