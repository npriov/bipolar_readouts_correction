#!/bin/bash

#example of running single bipolar acquisition correction function

# so the input of the function is the SVD you want to correct (5-7 dimensions),
#  a flag that takes values between 1-4 and indicates the preprocessing steps it will do
# (6 runs the whole pipeline)
# and a flag that indicates how to treat shots (none does no shot averaging). 
# For extra information, check the
# documentation 
input_svd="/raid/common/20181205_2/meas_MID212.complex.svd";
preprocessing_flag=4;
shots_flag="none";

matlab -r "bipolar_gradient_correction_one_acq('$input_svd','$preprocessing_flag','$shots_flag')"


# example of running two opposite direction bipolar acquisitions combination function

# so the input of the function is the SVDs you want to correct (5-7 dimensions),
#  a flag that takes values between 1-2 and indicates the preprocessing steps it will do
# (4 runs the whole pipeline)
# and a flag that indicates how to treat shots (none does no shot averaging). 
# For extra information, check the
# documentation 

input_svd1="/raid/common/20181205_2/meas_MID212.complex.svd";
input_svd2="/raid/common/20181205_2/meas_MID213.complex.svd";
preprocessing_flag=2;
shots_flag="none";

matlab -r "bipolar_gradient_correction_one_acq('$input_svd1','$input_svd2','$preprocessing_flag','$shots_flag')"
