# pulsarAnalysis #
Analysis scripts for crab pulsar data

Meant to be used on outputted files from Marten van Kerkwijk's scintellometry package:

https://github.com/mhvk/scintellometry

## GPs ##
Assorted tools for giant pulse analysis

### pulseFinder.py: ###

Finds giant pulses in a given time series of data. Stacks all input files.
Run as:

python pulseFinder.py foldspec1 foldspec2 ...

### pulseSpec.py: ###

Creates a dynamic spectrum for the largest giant pulse in a given time series of data. Stacks all input files.
Run as:

python pulseSpec.py foldspec1 foldspec2 ...

### dualPulseSpec.py: ###
Creates a side by side comparison of the dynamic spectra for the largest giant pulses in each of the two given time 
series of data.
Run as:

python dualPulseSpec.py foldspec1 foldspec2

### pulseProjFreq.py: ###

Projects the dynamic spectrum for the largest giant pulse in a given time series of data onto the frequency axis, and performs analysis of this spectrum. Stacks all input files.
Run as:

python pulseProjFreq.py foldspec1 foldspec2 ...

### pulseProjTime.py: ###

Projects the dynamic spectrum for the largest giant pulse in a given time series of data onto the time axis, and performs analysis of this profile. Stacks all input files.
Run as:

python pulseProjTime.py foldspec1 foldspec2 ...

### pulseCorr.py: ###
Calculates and plot the correlation coefficient between the spectra for each pair of pulses in all the files in the given directory.
Run as:

python pulseCorr.py /path/to/MP/files/ /path/to/IP/files/


## Misc ##
Assorted helper tools

### GMRTNaming.py ###
Converts between node/raw voltage numbers and dish names for the various GMRT antennae.
Run as:

python GMRTNaming.py dishName
or
python GMRTNaming.py node rawVolt

### GMRTStatus.py ###
Gets the current status (good/broken/missing) of the various GMRT antennae.
Run as:

python GMRTNaming.py dishName
or
python GMRTNaming.py node rawVolt

### GMRTDelay.py ###
Gets the delay in integer number of samples for the giant pulse at ~15:01:45 in the various GMRT antennae.
Run as:

python GMRTDelay.py dishName
or
python GMRTDelay.py node rawVolt

### voltToInt.py ###
Adds all input voltageDump files with associated phases together, and outputs the intensity associated with these files.
Run as:

python voltToInt.py foldspec1 foldspec2 ...
or
python voltToInt.py path/to/foldspecs/
