Eddy currents can result in phase and magnitude disturbances in MRI. While the additional phase is mostly constant in unipolar acquisitions, this can create a quite big phase difference between positive and negative echoes, if one is using a bipolar setup for speed gains/extra sampling. Here is for example the phase of a bipolar readout after unwrapping and regression of the phase evolution over the echo train. You can see clearly the linear phase accumulation:

![alt text](https://github.com/npriov/bipolar_readouts_correction/blob/master/raw_phase.png?raw=true)

This is easy to take out with a simple regression within slice. Turns out that the effect is also fairly stable over the whole echo train.

![alt text](https://github.com/npriov/bipolar_readouts_correction/blob/master/regr_phase.png?raw=true)

However, there are still a DC phase offset left for each phase direction (the stacked echo train looks quite jagged)

![alt text](https://github.com/npriov/bipolar_readouts_correction/blob/master/raw_phasete.png?raw=true)

This can also be corrected for quite easily with simple filtering. Note that, since this application was focused on T2* myelin mapping, hence sensitive to short T2*, we weighted the filter with a sigmoid funcion to the end of the echo train.

![alt text](https://github.com/npriov/bipolar_readouts_correction/blob/master/prcoc_phase.png?raw=true)

Some disturbances can still be seen, partially relating to the navigator and other scanner and sequency specifics; they can still be corrected with some spatial filtering. 
Note that magnitude also shows some disturbance, but much smaller.

A simple alternative to these corrections would be to measure with reversed gradient polarity and average. This also works pretty well, albeit at the cost of longer scans and assuming few phase differences between scans. 

Functions for both cases are implemented in the repo.

