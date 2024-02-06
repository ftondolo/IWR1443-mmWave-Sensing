import csv 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.signal import find_peaks
import numpy
import os

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Open CSV
with open('./output.csv') as file_obj: 
      
    # Create reader object
    reader_obj = csv.reader(file_obj) 

    #Variables to store data type  
    saved_state = []
    range_noise = 0
    range = 0
    noise = 0
    doppler = 0
    rangeArray = numpy.array([])
    dopplerArray = numpy.array([])
    for row in reader_obj: 
        # Classify input data based on first row ID tag
        if row[0]=="range_noise":
            range_noise = 1
            continue
        elif row[0]=="range":
            range = 1
            continue
        elif row[0]=="noise":
            noise = 1
            continue
        elif row[0]=="doppler":
            doppler = 1
            continue
        # Transform string to float array
        row = [float(num) for num in row]
        # If range and noise are being transmitted
        if range_noise:
            # Save state if none saved for next-iteration comparison
            if saved_state == []:
                saved_state = row
                continue
            plt.title("Range/Noise Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            # Adjust axis markings
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            plt.plot(row, color='red')
            # Identify and Plot maxima, height chosen semi-arbitrarily
            peaks, _ = find_peaks(row, height=9500)
            plt.plot(peaks, [row[i] for i in peaks], "x")
            plt.plot(saved_state, color='blue')
            plt.pause(0.2)
            plt.clf()
        # If range alone being tracked
        if range:
            plt.title("Range Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            # Adjust axis markings
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            # Identify and Plot maxima, height chosen semi-arbitrarily
            peaks, _ = find_peaks(row, height=9500)
            cs = plt.plot(row, color='blue')
            plt.plot(peaks, [row[i] for i in peaks], "x")
            plt.pause(0.2)
            plt.clf()
        # If noise alone being tracked
        if noise:
            plt.title("Noise Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            # Adjust axis markings
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            cs = plt.plot(row, color='red')
            plt.pause(0.2)
            plt.clf()
        # If doppler heatmap being received
        if doppler:
            # First save Range Array
            if len(rangeArray) == 0:
                rangeArray = numpy.array(row)
            # If Range Array saved then save Doppler Array
            elif len(dopplerArray) == 0:
                dopplerArray = numpy.array(row)
            # If Range and Doppler arrays have been saved, then received Range-Doppler
            else:
                # Reshape previously flattened array
                row = numpy.reshape(row, (16,256))
                fig = plt.figure()
                plt.title("Doppler Range Heatmap")
                plt.xlabel("Range (m)")
                plt.ylabel("Doppler (m/s)")
                cs = plt.contourf(rangeArray,dopplerArray,numpy.reshape(row, (16,256)))
                fig.colorbar(cs, shrink=0.9)
                plt.pause(0.2)
                plt.clf()
                plt.close()
                # Reset saved states for new iteration
                rangeArray = []
                dopplerArray = []