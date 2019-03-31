
%example of running single bipolar acquisition correction function
inp_list={'/raid/common/20181205_2/meas_MID212.complex.svd',...
    '/raid/common/20181108_1/meas_MID766.svd'};


for i=1:length(inp_list)
    % so the input of the function is the SVD you want to correct (5-7
    % dimensions),
    % a flag that takes values between 1-6 and indicates the preprocessing steps it will do. 
    % (6 runs the whole pipeline)
    % and a flag that tells it whether to average all shots before
    % correction ("all"), or correct some shots (for example first to second shot "1:2") 
    % or do the correction per shot ("none")
    % For extra information, check the documentation 
    ff=bipolar_gradient_correction_linearramp_filter(inp_list{i},"4","none");
end



%example of running two opposite direction bipolar acquisitions combination function

inp_list={{'/raid/common/20181205_2/meas_MID212.complex.svd','/raid/common/20181205_2/meas_MID213.complex.svd'},...
            {'/raid/common/20181108_1/meas_MID766.svd','/raid/common/20181108_1/meas_MID767.svd'}};

for i=1:length(inp_list)
    % so the input of the function are the SVD you want to correct (5-7
    % dimensions),
    % a flag that takes values between 1-2 and indicates the preprocessing steps it will do. 
    % (2 runs the whole pipeline)
    % and a flag that tells it whether to average all shots before
    % correction ("all"), or correct some shots (for example first to second shot "1:2") 
    % or do the correction per shot ("none")
    % For extra information, check the documentation 
    ff=opposite_bipolar_gradient_correction(inp_list{i}{1},inp_list{i}{2},"2","none");
end
