# Photoacoustic imaging with Pump Laser
This folder houses all the functions need to perform photoacoustic imaging using a pump laser and perform a streamlined analysis approach. This guide offers step-by-step instructions that a new Verasonics user can follow to perform PAI using a pump laser. The guide assumes the user has a lisenced version of verasonics installed and activated (if not, please follow the Verasonics Tutorial/Manuals to do so). The guide also assumes that the trigger connections between the Vantage and the pump laser is performed properly.    

## Setup and Imaging
This section guides the user to save the matlab variables required by the Verasonics software to perform PAI. The setup file is written for the L7-4 linear array transducer, which the user must change if a different transducer must be used. The section assumes that all the files in this repository are saved to the user's local system and added to current MATLAB path. The user can then follow the steps below to perform PAI, 

* Run "SetUpL7_4PA_Delay_BModeToggle.m" to save .mat file for VSX processing. A .mat file named "L7-4PA_PumpLaser.mat" will be saved under the VSX folder's MatFiles directory. 
```
>>> SetUpL7_4PA_Delay_BModeToggle
```
* Run saved .mat file using the Verasonics function "VSX"
```
>>> VSX; 
>>> Name of .mat file to process: L7-4PA_PumpLaser.mat
```
* Acquire the required PA data by pressing "Save Burst". "Next Bust" can be used to save a new set of bursts. By default, the acquired data will be saved under "PA_Acqs" directory under the Vantage software folder. Below is an example of a B-mode and PA image of a rat kidney. 
![An image of the VSX Gui with a sample B-mode & PA image of Kidney](https://github.com/VinVincely/BFIL_Verasonics_fns/blob/main/PA_PumpLaser/images/VSX_GUI.png)

* Ensure that each burst folder has a "info.txt" file. This will contain information of the current image acquistions such as wavelength, surface fluence, and other information that may be relavent (See image below for example).

* Once all the required bursts are collected, close the VSX Control GUI to stop the Vantage acquisition and return control to Matlab. 
![An image of the data acquisition folder](https://github.com/VinVincely/BFIL_Verasonics_fns/blob/main/PA_PumpLaser/images/burst_fldr.jpg)