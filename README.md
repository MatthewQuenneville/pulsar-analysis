# pulsar-analysis
Analysis scripts for crab pulsar data

Meant to be used on outputted files from Marten van Kerkwijk's scintellometry package:

https://github.com/mhvk/scintellometry

### pulseFinder.py: ###

Finds giant pulses in a given time series of data.
Run as:

python pulseFinder.py foldspec

### pulseSpec.py: ###

Creates a dynamic spectrum for the largest giant pulse in a given time series of data.
Run as:

python pulseSpec.py foldspec

### dualPulseSpec.py: ###
Creates a side by side comparison of the dynamic spectra for the largest giant pulses in each of the two given time 
series of data.
Run as:

python dualPulseSpec.py foldspec1 foldspec2

### pulseProjFreq.py: ###

Projects the dynamic spectrum for the largest giant pulse in a given time series of data onto the frequency axis, and performs analysis of this spectrum.
Run as:

python pulseProjFreq.py foldspec

### pulseProjTime.py: ###

Projects the dynamic spectrum for the largest giant pulse in a given time series of data onto the time axis, and performs analysis of this profile.
Run as:

python pulseProjTime.py foldspec

### pulseCorr.py: ###
Calculates and plot the correlation coefficient between the spectra for each pair of pulses in all the files in the given directory.
Run as:

python pulseCorr.py /path/to/MP/files/ /path/to/IP/files/
