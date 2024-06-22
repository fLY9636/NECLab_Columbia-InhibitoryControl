function nearestIDX = CrossSampling(set_loc, ref)
% CrossSampling search the index of an element in set_loc which is closest
% to the ref

[val, loc] = min(abs(set_loc-ref));
nearestIDX = loc(1);
end