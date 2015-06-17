# pulsar-analysis
Analysis scripts for crab pulsar data

Meant to be used on outputted files from Marten van Kerkwijk's scintellometry package:

https://github.com/mhvk/scintellometry

pulseFinder.py:
Finds giant pulses in a given time series of data.
Run as:

python pulseFinder.py foldspec

pulseSpec.py:
Creates a dynamic spectrum for the largest giant pulse in a given time series of data.
Run as:

python pulseSpec.py foldspec

dualPulseSpec.py:
Creates a side by side comparison of the dynamic spectra for the largest giant pulses in each of the two given time 
series of data.
Run as:

python dualPulseSpec.py foldspec1 foldspec2
