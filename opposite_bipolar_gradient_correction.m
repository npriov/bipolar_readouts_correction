function f= opposite_bipolar_gradient_correction(input_svd1,input_svd2,proc_steps, shots_flag)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this function receives two opposite polarity images
        
    % assumes that inputs are 5,6 or 7 dimensions SVD + a number between 1-2
    % indicating the preprocessing steps you would like to perform + a string 
    % to interpret how to treat shots. Possible options are "all" 
    % (average over shots), "none" (treats its shot independently), "num1:num2"
    % i.e. "1:2" which averages the first two shots 
    
    % the SVD is reformatted in (PE,FE,slices,shots,TE) before any
    % preoprocessing is done
    
    % corrected 5D SVDs will be created.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preprocessing steps 
    
    % 1. will just unwrap the phase, regress out increase along TE and then
    % do a complex combination as: 
    % sqrt(magnitude_pos^2+magnitude_neg^2)*exp(phase_pos+phase_neg)
    % per 10.1016/j.neuroimage.2018.11.040
    
    % 2. unwrap phase + low-pass filtering in in-plane. 
    % There still seems to be a 3D phase disturbance,
    % To preserve information in the TE direction, we apply a high pass
    % in plane filter (so orthogonal to echo plane)
    % The phase disturbance is not stable between TEs, so it is calculated
    % for each TE.
    % The precise threshold used can be set in the parameters.
 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %set up paths
    addpath('/misc/imeel/dezwart/matlab');
    addpath(genpath('/misc/imeel/priovoulosn2/matlab'));
    %addpath(genpath('~/Documents/MATLAB/bipolar_corrections'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    proc_steps=str2double(proc_steps);   
  
    if (proc_steps>2)
        f=0;
        return
    end

    
    %read in input svds
    data_pos=read_data(input_svd1);
    data_neg=read_data(input_svd2);   
    
    if (size(data_pos)~=size(data_neg))      
        error('size of the two input SVDs is not the same');
    end
    
    if (length(size(data))< 5 || length(size(data)) > 7)      
        error('size of the two input SVDs is not between 5 to 7 dimensions');
    end
        
    [a,b]=fileparts(input_svd1);
    b=strsplit(b,'.');
    
    if isfile(strcat(ls(strcat(a,'/',b{1},'*','prun/echo_times*'))))
        echo_times=read_data(strcat(ls(strcat(a,'/',b{1},'*','prun/echo_times*'))));
        echo_times=[echo_times(1).echo0.data',echo_times(1).echo1.data'];
    else
        error('cannot find echo time SVD. Check echo_times initialization');
    end
    % save_data(strcat(out_roi,'/echo_times.svd'),echo_times');
 
    %check data dimension and bring to 5D format (PE,FE,slices,shots,TE)
    if (length(size(data_pos))==6)
        data_pos=permute(squeeze(data_pos),[1,2,4,5,3]);  
        data_neg=permute(squeeze(data_neg),[1,2,4,5,3]);        
    end
    
    if (length(size(data_pos))==7)
        data_pos=mean(data_pos,7);
        data_pos=squeeze(data_pos);         
        data_neg=mean(data_neg,7);
        data_neg=squeeze(data_neg);  
    end  

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parameter setup
   
    %set up names for output file of each preprocessing step.
    [a,b,~]=fileparts(input_svd1);  
    out_comb=strcat(a,"/",b,"_comb.svd");
    out_comb_2dfilt=strcat(a,"/",b,"_comb_2dfilt.svd");
   
    
    %define highpass for filter
    spatial_threshold_freq=40; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
  
    % unwrap phase
    unwrap_phase_pos=unwrap(angle(data_pos),[],5);
    phase_increase_pos=zeros(size(unwrap_phase_pos));
    unwrap_phase_neg=unwrap(angle(data_neg),[],5);
    phase_increase_neg=zeros(size(unwrap_phase_neg));

    % regress out linear increase along TE
    for i = 1:size(unwrap_phase_pos,1)
        for j = 1:size(unwrap_phase_pos,2)
            for k = 1:size(unwrap_phase_pos,3)
                for u = 1:size(unwrap_phase_pos,4)
                    p=polyfit(echo_times',squeeze(unwrap_phase_pos(i,j,k,u,:)),1);
                    yfit=polyval(p, echo_times');
                    unwrap_phase_pos(i,j,k,u,:)=squeeze(unwrap_phase_pos(i,j,k,u,:))-yfit;
                    phase_increase_pos(i,j,k,u,:)=yfit;
                    
                    p=polyfit(echo_times',squeeze(unwrap_phase_neg(i,j,k,u,:)),1);
                    yfit=polyval(p, echo_times');
                    unwrap_phase_neg(i,j,k,u,:)=squeeze(unwrap_phase_neg(i,j,k,u,:))-yfit;
                    phase_increase_neg(i,j,k,u,:)=yfit;                    
                end
            end

        end
    end

    %complex combination of polarities
    
    %interpret shot average flag
    if (shots_flag=="all")                
        unwrap_phase_pos=mean(unwrap_phase_pos,4);
        unwrap_phase_neg=mean(unwrap_phase_neg,4);
        mag_both_acq=sqrt(mean(abs(data_pos),4).^2+mean(abs(data_neg),4).^2);

    elseif (shots_flag=="none")        
        b=1;
    else            
        unwrap_phase_pos=mean(unwrap_phase_pos(:,:,:,eval(shots_flag),:),4);
        unwrap_phase_neg=mean(unwrap_phase_neg(:,:,:,eval(shots_flag),:),4);
        mag_both_acq=sqrt(mean(abs(data_pos(:,:,:,eval(shots_flag),:)),4).^2+mean(abs(data_neg(:,:,:,eval(shots_flag),:)),4).^2);

    end
  
    mag_both_acq_normed=mag_both_acq./mag_both_acq(:,:,:,:,1);
    mag_both_acq_normed(isnan(mag_both_acq_normed))=0;   
    phase_both_acq=unwrap_phase_pos+unwrap_phase_neg;
    
    %save output and exit if needed
    save_data(out_comb,mag_both_acq_normed.*exp((1i)*(phase_both_acq))); 
    if (proc_steps==1)
        f=1;
        return
    end

    % spatial filter for each TE and slice
    phase_both_acq_ramped_filt_spatialfilt=zeros(size(phase_both_acq));
    for (slice=1:size(phase_both_acq_ramped_filt_spatialfilt,3))
        for (shot=1:size(phase_both_acq_ramped_filt_spatialfilt,4))
            for (echo=1:size(phase_both_acq_ramped_filt_spatialfilt,5))
                phase_both_acq_ramped_filt_spatialfilt(:,:,slice,shot,echo)=squeeze(phase_both_acq(:,:,slice,shot,echo))...
                +abs((ifft(bhp((fft(squeeze(phase_both_acq(:,:,slice,shot,echo)))),spatial_threshold_freq,3))));
            end
        end
    end
    
    save_data(out_comb_2dfilt,mag_both_acq_normed.*exp((1i)*(phase_both_acq_ramped_filt_spatialfilt)));
    f=1;
    return
   
    
end












