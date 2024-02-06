# IWR1443-mmWave-Sensing
A Python script to graph and save data from a IWR1443Boost dev board 

## How It Works
1. Plug in your IWR1443 Development Board, this repository specifically has been tested with the IWR143Boost board running SDK 2.1
2. Specify your UART and Data ports in the code `serialConfig` function
3. Run the code using `python3 data_plotting.py`

## Changing Display Type
The code automatically detects whether the IWR1443 is attempting to display a range profile, a noise profile, both a range and noise profile, an azimuth heatmap, a doppler ehatmap, or statistics, it is however limited to reading only one of these options at any given time. To specify, you must chaneg your profile.cfg file. Specifically, the `guiMonitor` line where each subsequent digit 0/1 delineates whether that particular GUI function is activated or deactivated. Consult the code, particularly the function `parseConfig`, to determine which of these to flip to display the desired data. 

## Visualizing Saved Data
The data being visualized if it is one of the supported functions (range, noise, range/noise, doppler) will be saved to an `output.csv` file in the code's directory which can be later read with `readFile.py`

