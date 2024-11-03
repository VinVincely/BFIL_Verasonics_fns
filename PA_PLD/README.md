# Photoacoustic imaging with Pump Laser
This folder houses all the functions needed to perform photoacoustic imaging using a pulsed laser diode array triggered by a laser driver. 

## Setup and Imaging
This section guides the user to save the matlab variables required by the Verasonics software using the associated acquisiton scripts. The acquisition script is set up to perform high frame rate imaging (up to ~5000 Hz) where real time image visualization looks terribly "jerky". To reconstruct all data, the associated "RECON" script must be utilized. The setup files are written for the L7-4 linear array transducer, which the user must change if a different transducer must be used. The user can then follow the steps below to perform PAI, 

***** NOTE: UNDER CONSTRUCTION! *****

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

## Loading and Visualization of acquired data
This section guides the user to load the acquired data and saving them to easily accessible variables. The loaded data can then be easily visualized using the Visualization functions in this repository. 

* Load acquired verasonics data using "load_Verasonics_Bursts". NOTE: The present working directory must contain all the acquired bursts.
```
>>> [PA_data, Bmode_data, burst_info] = load_Verasonics_Bursts('.');
```
Here, "PA_data" & "Bmode_data" are cells that contain the acquired I/Q data of the photoacoustic and ultrasound B-mode images (Size of cell will be equal to number of bursts). "burst_info" stores all the information cataloged in the "info.txt" files within each burst folder. 

* To convert the I/Q data to image files, the function "get_image_from_VSX_data" can be used as follows, 
```
>>> PA_data = get_image_from_VSX_data(PA_data, 1);
```

* If the user needs to retrieve specific frames from the entire PA/US acquisition ensemble, the function "get_nAcq_VSX_Bursts" as follows, 
```
>>> PA_data = get_nAcq_VSX_Bursts(PA_data, [1 5 9 13 17]);   
```
Here, specific frames numbered - 1, 5, 9, 13 & 17 are retrieved for each burst and saved to the "PA_data" cell. 

* To visualize the data, the user can use the function "visualize_VSX_sPAI_data". The function can be used to visualize multiple bursts simultaneously. NOTE: Each bust must only contain one frame for code to function.  
```
visualize_VSX_sPAI_data(PA_data, Bmode_data, meanFlnc, lambdas, [1 14]); 
```
The above function will plot the 1st and 14th burst (as seen below). 
![An image of the data acquisition folder](https://github.com/VinVincely/BFIL_Verasonics_fns/blob/main/PA_PumpLaser/images/visualizePA.png)