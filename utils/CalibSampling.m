function idx = CalibSampling(PTime, Pupil_in_imec, t)
% CalibSampling takes the inputs PTime (timestamp of triggers sent by XPC),
% Pupil_in_imec (XPC triggers' indices in the imec timestamp) and t (time of interest in the behavior data). It outputs idx, the index
% of t in the imec timestamp.
if length(t) == 1
    Pidx = CrossSampling(PTime, t);                                   % index in the pupil trigger timestamp
    dT = t-PTime(Pidx);                                                         % time difference (s) between the t and nearest pupil trigger
    idx = Pupil_in_imec(Pidx)+round(dT*30000);                           % get the corresponding index of that pupil trigger in the imec time and add the time difference (a.p. fs converted)
else
    idx = [];
    for i = 1:length(t)
        Pidx = CrossSampling(PTime, t(i));
        dT = t(i)-PTime(Pidx);
        idx(end+1, 1) = Pupil_in_imec(Pidx)+round(dT*30000);
    end
end
end
