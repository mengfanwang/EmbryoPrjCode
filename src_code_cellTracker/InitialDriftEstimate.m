function drift = InitialDriftEstimate(time_ordered_data)
% estimate the drift by aligning the data with translation only alignemnt.

ndim = ndims(time_ordered_data{1});
drift = zeros(numel(time_ordered_data), ndim);
[optimizer, metric] = imregconfig('Monomodal');
for i=2:numel(time_ordered_data)
    disp(i);
    tform = imregtform(time_ordered_data{i}, time_ordered_data{i-1},...
        'translation',optimizer,metric);
    drift(i,:) = -tform.T(ndim+1, 1:ndim);
end
end