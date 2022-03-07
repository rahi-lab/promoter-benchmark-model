function average = calculate_average(means)

N_columns = size(means, 2);
N_rows = size(means, 1);

average = [];

for i = 1:N_columns
    sum = 0;
    count = 0;
    %means might contain NaN values so we have to check for that
    %Matlab's mean(array) returns NaN if there are NaN values in array
    for j = 1:N_rows
        if(~isnan(means(j,i)))
                sum = sum + means(j,i);
                count = count + 1;
        end
    end
    average = [average, sum/count];
end

    