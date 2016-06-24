function [t_corr,t_corr_t] = correlation_time(data, t_max)
%t_max can't be greater than the rows of data

size_data = size(data);

%compute variance of data
var_data = var(data);
varDatat = repmat(var_data,t_max,1);

%compute square expectation value
data_expected = mean(data);
data_expected_sq = data_expected.^2;
data_expected_sq = repmat(data_expected_sq,t_max,1);

%compute correlation expectation value
corr_sq_expected = zeros(t_max,size_data(2));
for t = 1:t_max
    if(mod(t,100) == 0)
        t
    end
    corr_sq_sum = 0;
    for i = 1:(size_data(1) - t)
        corr_sq_sum = corr_sq_sum + data(i,:).*data(i + t,:);
    end

    corr_sq_expected(t,:) = corr_sq_sum/(size_data(1) - t);
end

t_corr_t = corr_sq_expected - data_expected_sq;
t_corr_t = t_corr_t./varDatat;
t_corr = sum(corr_sq_expected - data_expected_sq);
t_corr = t_corr./var_data;
