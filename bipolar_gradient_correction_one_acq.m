function f= bipolar_gradient_correction_one_acq(input_svd,proc_steps, shots_flag)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this function receives as input bipolar readouts and tries to reduce
    % the eddy current effect on the phase. 
    % assumes that input is 5,6 or 7 dimensions SVD + a number between 1-4
    % indicating the preprocessing steps you would like to perform + a string 
    % to interpret how to treat shots. Possible options are "all" 
    % (average over shots), "none" (treats its shot independently), "num1:num2"
    % i.e. "1:2" which averages the first two shots 
    
    % the SVD is reformatted in (PE,FE,slices,shots,TE) before any
    % preoprocessing is done
    
    % corrected 5D SVDs will be created.
    % In the output, the magnitude of each echo train is normalized with the magnitude of
    % the first echo.
    
    % the function will save a file for each preprocessing step it gets to
    % do, e.g. if you input "3", it will output all 3 steps. This is done
    % to facilitate comparisons.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % preprocessing steps 
    
    % 1. raw output, will just unwrap the phase and regress out increase along TE. 
    % I put it there tofaciliate comparisons
    
    % 2. unwrap phase + linear ramp. This is based on the algorithm originally 
    % described by Lu et al.,2008 https://doi.org/10.1002/mrm.21583 to
    % reduce the dominant eddy current effect.
    % The implementation here calculates the ramp for each TE and each
    % slice for a few lines in the center and then uses the median across 
    % TEs with the same polarity; 
    % The ramp seems to be pretty stable across TE, though there is some 
    % systematic difference; in practice results are more robust if the
    % median ramp is used so that's the default here. This behavior can be
    % changed in the parameter setup.
    % The ramp is evident only within brain, so we look for it only in a
    % subsample of the data (defined as acq_space,pha_space)
        % To calculate the ramp, we average across a few lines to get a robust signal
    % and then fit a linear regressor
    
    % 3. unwrap phase + linear ramp + band-pass filtering in TE. A sizeable phase
    % disturbance due to eddy currents is still evident that has the form
    % of a saw tooth oscillator, with a period of 2*(time between TEs); we 
    % therefore filter out high frequencies with a non-causal filter 
    % (idealfilter function from matlab).
    % To account for the navigator, which makes sampling irregular, the
    % data are linearly interpolated.
    % Assuming that the fast decay of myelin is of interest, the filter is not
    % applied in the initial TEs; the data are combined with a sigmoid.
    
    % 4. unwrap phase + linear ramp + band-pass filtering in TE + high-pass
    % filtering in space. 
    % Eddy currents seem to leave also a 3D phase disturbance 
    % (Yu et al., 2010,  https://doi.org/10.1002/jmri.22111)
    % To preserve information in the echo-train direction, we apply a high pass
    % filter orthogonal to echo-train
    % The phase disturbance is not stable between TEs, so it is calculated
    % for each TE.
    % The precise threshold used can be set in the parameters.
    
    % 5. unwrap phase + linear ramp + band-pass phase filtering in TE 
    % + band-pass magnitude filtering in TE
    % Magnitude also shows a sawtooth pattern. Here, we regress out a
    % single exponential decay in the magnitude and apply a filter for the 
    % magnitude, exactly the same as we do in step 3.
    % We then combine with the phase output of step 3.
    
    % 6. unwrap phase + linear ramp + band-pass phase filtering in TE 
    % + band-pass magnitude filtering in TE + high-pass phase filtering 
    % in space. 
    % Similar to step 4, but we use the band-pass magnitude filtered data. 
    % instead of the non-bandpass filtered data.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %set up paths
    addpath('/misc/imeel/dezwart/matlab');
    addpath(genpath('/misc/imeel/priovoulosn2/matlab'));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    proc_steps=str2double(proc_steps);   
 
    if (proc_steps>6)
        f=0;
        return
    end


    %read in input svd
    data=read_data(input_svd);
    [a,b]=fileparts(input_svd);
    b=strsplit(b,'.');
    if isfile(strcat(ls(strcat(a,'/',b{1},'*','prun/echo_times*'))))
        echo_times=read_data(strcat(ls(strcat(a,'/',b{1},'*','prun/echo_times*'))));
        echo_times=[echo_times(1).echo0.data',echo_times(1).echo1.data'];
    else
        error('cannot find echo time SVD. Check echo_times initialization');
    end
    % save_data(strcat(out_roi,'/echo_times.svd'),echo_times');
    
    if (length(size(data))< 5 || length(size(data)) > 7)      
        error('size of the input SVD is not between 5 to 7 dimensions');
    end

    %check data dimension and bring to 5D format (PE,FE,slices,shots,TE)
    if (length(size(data))==6)
        data=permute(squeeze(data),[1,2,4,5,3]);  
    end
    
    if (length(size(data))==7)
        data=mean(data,7);
        data=squeeze(data);  
    end  
   
    %define a rough brain mask based on amplitude
    rr=abs(data);
    rr=(rr(:,:,:,1,3)-rr(:,:,:,1,size(echo_times,2)/2));
    rr(rr<20)=0;
    rr(rr>0)=1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % parameter setup
   
    %set up names for output file of each preprocessing step.
    [a,b,~]=fileparts(input_svd);
    out_raw=strcat(a,"/",b,"_raw.svd");
    out_ramp=strcat(a,"/",b,"_ramp.svd");
    out_ramp_filt=strcat(a,"/",b,"_rampedfilt.svd");
    out_ramp_filt_3d=strcat(a,"/",b,"_rampedfilt3d.svd");
    out_ramp_filt_mag_filt=strcat(a,"/",b,"_rampedfilt_magfilt.svd");   
    out_ramp_filt_mag_filt_3d=strcat(a,"/",b,"_rampedfilt_magfilt_3d.svd");
    

    %linear ramp parameters
    % data space within which we will look for a linear ramp; 
    % needs to be within what we are imaging so we can detect the ramp.
    pha_space=ceil(size(data,2)/2)-5:ceil(size(data,2)/2)+5;
    
    temp=(rr(:,pha_space,:));
    %find minimum extension of brain mask within that pha_space
    brain_extent=floor((min(min(sum(temp>0,1)))-5)/2);
    clear temp
    clear rr
    %acq_space=(ceil(size(data,1)/2)-15):(ceil(size(data,1)/2)+15);
    acq_space=(ceil(size(data,1)/2)-brain_extent):(ceil(size(data,1)/2)+brain_extent);

    
    % this flag controls if we will average the calculated ramp across TEs
    % of the same polarity. In practice, while there is some structured
    % ramp variation between TEs, so there is an argument for putting this
    % flag to 0, it seemed that in some scans that include a lot of phase
    % instability, i.e. CSF, some fits for some TEs were bad. If one wants
    % to inspect the ramps, you can plot the variable angle_evolution_with_te
    average_ramp_across_TEs=1;
    
    %filtering parameters
    % echo from which on it is assumed that the myelin signal of the T2*
    % decay has relaxed completely (and therefore after this, we apply 
    % the high pass filter )
    sigm_step=17;
    
    %define bandpass for filter
    high_freq=(1/(2*(echo_times(2)-echo_times(1)))-0.1); %this is kHz, assuming input echo times are in ms
    interval = [0 high_freq]; 
    
    %define highpass for spatial filter
    spatial_threshold_freq=40; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %get normalized magnitude
    abs_normed=abs(data)./abs(data(:,:,:,:,1));
    abs_normed(isnan(abs_normed))=0;
        
    % unwrap phase
    unwrap_phase=unwrap(angle(data),[],5);
    phase_increase=zeros(size(unwrap_phase));

    % regress out linear increase along TE
    for i = 1:size(unwrap_phase,1)
        for j = 1:size(unwrap_phase,2)
            for k = 1:size(unwrap_phase,3)
                for u = 1:size(unwrap_phase,4)
                    p=polyfit(echo_times',squeeze(unwrap_phase(i,j,k,u,:)),1);
                    yfit=polyval(p, echo_times');
                    unwrap_phase(i,j,k,u,:)=squeeze(unwrap_phase(i,j,k,u,:))-yfit;
                    phase_increase(i,j,k,u,:)=yfit;
                end
            end

        end
    end
    
    % regress out exponential decay of amplitude
    abs_neg=abs(data);
    resid_neg=zeros(size(data));
    main=zeros(size(data));
    for i = 1:size(abs_neg,1)
        for j = 1:size(abs_neg,2)
            for k = 1:size(abs_neg,3)
                for u = 1:size(abs_neg,4)

                    p=polyfit(echo_times',log(squeeze(abs_neg(i,j,k,u,:))),1);
                    yfit=polyval(p, echo_times');
                    resid_neg(i,j,k,u,:)=exp(log(squeeze(abs_neg(i,j,k,u,:)))-yfit);
                    main(i,j,k,u,:)=exp(log(squeeze(abs_neg(i,j,k,u,:))));
                end
            end

        end
    end
    resid_neg(isnan(resid_neg))=0;

    %interpret shot average flag
    if (shots_flag=="all")        
        abs_normed=(mean(abs_normed,4));
        unwrap_phase=mean(unwrap_phase,4);
        resid_neg=mean(resid_neg,4);
        main=mean(main,4);
    elseif (shots_flag=="none")        
        b=1;
    else    
        abs_normed=mean(abs_normed(:,:,:,eval(shots_flag),:),4);
        unwrap_phase=mean(unwrap_phase(:,:,:,eval(shots_flag),:),4);
        resid_neg=mean(resid_neg(:,:,:,eval(shots_flag),:),4);
        main=mean(main(:,:,:,eval(shots_flag),:),4);
    end
    
    %save output and exit if needed
    save_data(out_raw,abs_normed.*exp((1i)*(unwrap_phase)));
    if (proc_steps==1)
        f=1;
        return
    end
    
      
    %apply processing
    %linear phase ramp
    angle_evolution_with_te=zeros(size(unwrap_phase,5),1);
    unwrap_phase_rampcor=zeros(size(unwrap_phase));
    
    % calculate linear ramp fit for each slice, TE and average
    for k = 1:size(unwrap_phase,3)
        for u = 1:size(unwrap_phase,4)
            for f= 1:size(unwrap_phase,5)                               
                pos=mean(squeeze(unwrap_phase(acq_space,pha_space,k,u,f)),2);
                p1 = robustfit(1:length(pos),pos');  
                angle_evolution_with_te(f)=p1(2);
            end
            
            %check if you should use the median of the ramp fits
            if (average_ramp_across_TEs==1)
                angle_evolution_with_te(angle_evolution_with_te>0)=median(angle_evolution_with_te(angle_evolution_with_te>0));
                angle_evolution_with_te(angle_evolution_with_te<0)=median(angle_evolution_with_te(angle_evolution_with_te<0));
            end
            
            %correct the data          
            for f= 1:size(unwrap_phase,5)                                        
                tmp=(angle_evolution_with_te(f)*(1-size(unwrap_phase,1)/2:size(unwrap_phase,1)/2));
                tmp=repmat(tmp, size(unwrap_phase,2),1)';
                unwrap_phase_rampcor(:,:,k,u,f)=unwrap_phase(:,:,k,u,f)-tmp;   
            end
        end
    end

    
    %save output and exit if needed
    save_data(out_ramp,abs_normed.*exp((1i)*(unwrap_phase_rampcor)));
    if (proc_steps==2)
        f=1;
        return
    end
    
    %filter along TE
    % define a weighting function to avoid the initial myelin related change 
    weight_function=sigmoid(echo_times,echo_times(sigm_step))';
    weight_function_inverted=abs(weight_function-1);

    phase_one_shot_ramped1=permute(unwrap_phase_rampcor, [5 1 2 3 4]);
    phase_one_shot_ramped_end=phase_one_shot_ramped1.*weight_function;
    phase_one_shot_ramped_end=permute(phase_one_shot_ramped_end, [2 3 4 5 1]);
    phase_one_shot_ramped_filt=zeros(size(phase_one_shot_ramped_end));
    
    abs_normed1=permute(resid_neg, [5 1 2 3 4]);
    abs_normed_end=abs_normed1.*weight_function;
    abs_normed_end=permute(abs_normed_end, [2 3 4 5 1]);
    abs_normed_end_filt=zeros(size(abs_normed_end));

    %fitting for each slice, TE and average
    for k = 1:size(phase_one_shot_ramped_end,1)
        for u = 1:size(phase_one_shot_ramped_end,2)
            for f= 1:size(phase_one_shot_ramped_end,3)
                for h= 1:size(phase_one_shot_ramped_end,4)

                    y=timeseries(squeeze(phase_one_shot_ramped_end(k,u,f,h,:)),echo_times);
                    tempf=idealfilter(y,interval,'pass') + mean(y);
                    phase_one_shot_ramped_filt(k,u,f,h,:) = tempf.Data;

                    y=timeseries(squeeze(abs_normed_end(k,u,f,h,:)),echo_times);
                    tempf=idealfilter(y,interval,'pass') + mean(y);
                    abs_normed_end_filt(k,u,f,h,:) = tempf.Data;                   
                    
                end
            end
        end
    end
    
    phase_one_shot_ramped_start=permute(unwrap_phase_rampcor, [5 1 2 3 4]);
    phase_one_shot_ramped_start=phase_one_shot_ramped_start.*weight_function_inverted;
    phase_one_shot_ramped_start=permute(phase_one_shot_ramped_start,[2 3 4 5 1]);
    phase_one_shot_ramped_filt=phase_one_shot_ramped_filt+phase_one_shot_ramped_start;
    
    abs_normed_start=permute(resid_neg, [5 1 2 3 4]);
    abs_normed_start=abs_normed_start.*weight_function_inverted;
    abs_normed_start=permute(abs_normed_start,[2 3 4 5 1]);
    abs_normed_end_filt=abs_normed_end_filt+abs_normed_start;

    %save output and exit if needed
    save_data(out_ramp_filt,abs_normed.*exp((1i)*(phase_one_shot_ramped_filt)));
    if (proc_steps==3)
        f=1;
        return
    end
    
    % spatial filter for each TE and slice
    phase_one_shot_ramped_filt_spatialfilt=zeros(size(phase_one_shot_ramped_end));
    for (slice=1:size(phase_one_shot_ramped_filt,3))
        for (shot=1:size(phase_one_shot_ramped_filt,4))
            for (echo=1:size(phase_one_shot_ramped_filt,5))
                phase_one_shot_ramped_filt_spatialfilt(:,:,slice,shot,echo)=squeeze(phase_one_shot_ramped_filt(:,:,slice,shot,echo))+abs((ifft(bhp((fft(squeeze(phase_one_shot_ramped_filt(:,:,slice,shot,echo)))),spatial_threshold_freq,3))));
            end
        end
    end
    
    %save output and exit if needed
    save_data(out_ramp_filt_3d,abs_normed.*exp((1i)*(phase_one_shot_ramped_filt_spatialfilt)));
    if (proc_steps==4)
        f=1;
        return
    end

    % save bandpass filtered amplitude residuals along TE
    main=abs_normed_end_filt+main;
    main=abs(main)./abs(main(:,:,:,:,1));
    main(isnan(main))=0;
    
    save_data(out_ramp_filt_mag_filt,main.*exp((1i)*(phase_one_shot_ramped_filt)));
    if (proc_steps==5)
        f=1;
        return
    end    

    save_data(out_ramp_filt_mag_filt_3d,main.*exp((1i)*(phase_one_shot_ramped_filt_spatialfilt)));
    f=1;
    return    
       
end






