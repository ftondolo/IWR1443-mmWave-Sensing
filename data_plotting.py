import serial
import time
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtGui, QtCore
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.ticker as ticker
from operator import add
import numpy as np
import math
import os
from itertools import chain
from scipy.signal import find_peaks
import csv

abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)

# Change the configuration file name
configFileName = './profile.cfg'

CLIport = {}
Dataport = {}
byteBuffer = np.zeros(2**15,dtype = 'uint8')
byteBufferLength = 0
rangeAzimuthHeatMapGridInit = 0


# ------------------------------------------------------------------

# Function to configure the serial ports and send the data from
# the configuration file to the radar
def serialConfig(configFileName):
    """
    Configures the serial ports and sends the CLI commands to the radar. 
    It outputs the serial objects for the data and CLI ports.
    """
    
    global CLIport
    global Dataport
    # Open the serial ports for the configuration and the data ports
    
    # Linux
    #CLIport = serial.Serial('/dev/ttyACM0', 115200)
    #Dataport = serial.Serial('/dev/ttyACM1', 921600)

    # Windows
    #CLIport = serial.Serial('/dev/COM4', 115200)
    #Dataport = serial.Serial('/dev/COM5', 921600)
    
    # Mac
    CLIport = serial.Serial('/dev/tty.usbmodemR10310411', 115200)
    Dataport = serial.Serial('/dev/tty.usbmodemR10310414', 921600)

    # Read the configuration file and send it to the board
    config = [line.rstrip('\r\n') for line in open(configFileName)]
    for i in config:
        CLIport.write((i+'\n').encode())
        print(i)
        time.sleep(0.01)
        
    return CLIport, Dataport

# ------------------------------------------------------------------
# Global variables to store selected GUI representation
det_obj = 0
range_prof = 0
noise_prof = 0
azimuth_heat = 0
doppler_heat = 0
stats = 0
# Function to parse the data inside the configuration file
def parseConfigFile(configFileName):
    global det_obj, range_prof, noise_prof, azimuth_heat, doppler_heat, stats
    """
     Parses the configuration file to extract the configuration parameters. 
     It returns the configParameters dictionary with the extracted parameters.
    """
    configParameters = {} # Initialize an empty dictionary to store the configuration parameters
    
    # Read the configuration file and send it to the board
    config = [line.rstrip('\r\n') for line in open(configFileName)]
    for i in config:
        
        # Split the line
        splitWords = i.split(" ")
        
        # Hard code the number of antennas, change if other configuration is used
        numRxAnt = 4
        numTxAnt = 2
        
        # Get the information about the profile configuration
        if "profileCfg" in splitWords[0]:
            startFreq = int(float(splitWords[2]))
            idleTime = int(splitWords[3])
            rampEndTime = float(splitWords[5])
            freqSlopeConst = float(splitWords[8])
            numAdcSamples = int(splitWords[10])
            numAdcSamplesRoundTo2 = 1
            
            while numAdcSamples > numAdcSamplesRoundTo2:
                numAdcSamplesRoundTo2 = numAdcSamplesRoundTo2 * 2
            digOutSampleRate = int(splitWords[11])
            
        # Get the information about the frame configuration    
        elif "frameCfg" in splitWords[0]:
            chirpStartIdx = int(splitWords[1]);
            chirpEndIdx = int(splitWords[2]);
            numLoops = int(splitWords[3]);
            numFrames = int(splitWords[4]);
            framePeriodicity = int(splitWords[5]);
        elif "guiMonitor" in splitWords[0]:
            det_obj = int(splitWords[1])
            range_prof = int(splitWords[2])
            noise_prof = int(splitWords[3])
            azimuth_heat = int(splitWords[4])
            doppler_heat = int(splitWords[5])
            stats = int(splitWords[6])
            
    # Combine the read data to obtain the configuration parameters           
    numChirpsPerFrame = (chirpEndIdx - chirpStartIdx + 1) * numLoops
    configParameters["numDopplerBins"] = numChirpsPerFrame / numTxAnt
    configParameters["numRangeBins"] = numAdcSamplesRoundTo2
    configParameters["rangeResolutionMeters"] = (3e8 * digOutSampleRate * 1e3) / (2 * freqSlopeConst * 1e12 * numAdcSamples)
    configParameters["rangeIdxToMeters"] = (3e8 * digOutSampleRate * 1e3) / (2 * freqSlopeConst * 1e12 * configParameters["numRangeBins"])
    configParameters["dopplerResolutionMps"] = 3e8 / (2 * startFreq * 1e9 * (idleTime + rampEndTime) * 1e-6 * configParameters["numDopplerBins"] * numTxAnt)
    configParameters["maxRange"] = (300 * 0.9 * digOutSampleRate)/(2 * freqSlopeConst * 1e3)
    configParameters["maxVelocity"] = 3e8 / (4 * startFreq * 1e9 * (idleTime + rampEndTime) * 1e-6 * numTxAnt)
    
    return configParameters # returns dictionary containing config key values
   
# ------------------------------------------------------------------

# Funtion to read and parse the incoming data
def readAndParseData14xx(Dataport, configParameters):
    """
    It reads the data from the data serial port and parses the recived buffer to extract the data 
    from the Detected objects package only. Othe package types (range profile, range-azimuth heat map...) 
    could be done similarly but have not been implemented yet. This functions returns a boolean variable 
    (dataOK) that stores if the data has been correctly, the frame number and the detObj dictionary with 
    the number of detected objects, range (m), doppler velocity (m/s), peak value and 3D position (m).
    """
    global byteBuffer, byteBufferLength
    
    # Constants
    NUM_ANGLE_BINS = 64
    OBJ_STRUCT_SIZE_BYTES = 12
    BYTE_VEC_ACC_MAX_SIZE = 2**15
    MMWDEMO_UART_MSG_DETECTED_POINTS = 1
    MMWDEMO_UART_MSG_RANGE_PROFILE   = 2
    MMWDEMO_UART_MSG_NOISE_PROFILE   = 3
    MMWDEMO_UART_MSG_AZIMUT_STATIC_HEAT_MAP = 4
    MMWDEMO_UART_MSG_RANGE_DOPPLER   = 5
    MMWDEMO_UART_MSG_STATS = 6
    maxBufferSize = 2**15
    magicWord = [2, 1, 4, 3, 6, 5, 8, 7]
    
    # Initialize variables
    magicOK = 0 # Checks if magic number has been read
    dataOK = 0 # Checks if the data has been read correctly
    frameNumber = 0
    detObj = {}
    
    readBuffer = Dataport.read(Dataport.in_waiting) # read buffer into readBuffer from serial data port
    byteVec = np.frombuffer(readBuffer, dtype = 'uint8') # store read data from buffer in var byteVec
    byteCount = len(byteVec)
    
    # byteBufferLength = 0
    # Check that the buffer is not full, and then add the data to the buffer and update byteBufferLength if requiered
    if (byteBufferLength + byteCount) < maxBufferSize:
        byteBuffer[byteBufferLength:byteBufferLength + byteCount] = byteVec[:byteCount]
        byteBufferLength = byteBufferLength + byteCount
        
    # Check that the buffer has some data
    if byteBufferLength > 16: # apparently it has to have at least 16 bytes in buffer
        
        # Check for all possible locations of the magic word
        possibleLocs = np.where(byteBuffer == magicWord[0])[0] # It returns a tuple of indices if an only condition is given, the indices where the condition is True.

        # Confirm that is the beginning of the magic word and store the index in startIdx
        startIdx = []
        for loc in possibleLocs:
            check = byteBuffer[loc:loc+8]
            if np.all(check == magicWord):
                startIdx.append(loc)
               
        # Check that startIdx is not empty
        if startIdx:
            
            # Remove the data before the first start index
            if startIdx[0] > 0 and startIdx[0] < byteBufferLength:
                byteBuffer[:byteBufferLength-startIdx[0]] = byteBuffer[startIdx[0]:byteBufferLength]
                byteBuffer[byteBufferLength-startIdx[0]:] = np.zeros(len(byteBuffer[byteBufferLength-startIdx[0]:]),dtype = 'uint8')
                byteBufferLength = byteBufferLength - startIdx[0]
                
            # Check that there have no errors with the byte buffer length
            if byteBufferLength < 0:
                byteBufferLength = 0
                
            # word array to convert 4 bytes to a 32 bit number
            word = [1, 2**8, 2**16, 2**24]
            
            # Read the total packet length
            totalPacketLen = np.matmul(byteBuffer[12:12+4],word) # Matrix product of two arrays
            # total packet lenght stored in 12:12+4 position in 4 bytes in the packet
            
            # Check that all the packet has been read
            if (byteBufferLength >= totalPacketLen) and (byteBufferLength != 0):
                magicOK = 1
    
    # If magicOK is equal to 1 then process the message
    if magicOK:
        # word array to convert 4 bytes to a 32 bit number
        word = [1, 2**8, 2**16, 2**24]
        
        # Initialize the pointer index
        idX = 0
        
        # Read the header
        magicNumber = byteBuffer[idX:idX+8]
        idX += 8
        version = format(np.matmul(byteBuffer[idX:idX+4],word),'x')
        idX += 4
        totalPacketLen = np.matmul(byteBuffer[idX:idX+4],word)
        idX += 4
        platform = format(np.matmul(byteBuffer[idX:idX+4],word),'x')
        idX += 4
        frameNumber = np.matmul(byteBuffer[idX:idX+4],word)
        idX += 4
        timeCpuCycles = np.matmul(byteBuffer[idX:idX+4],word)
        idX += 4
        numDetectedObj = np.matmul(byteBuffer[idX:idX+4],word)
        idX += 4
        numTLVs = np.matmul(byteBuffer[idX:idX+4],word) # number of Data Structures in package (4 bytes)
        idX += 4
        
        # Read the TLV messages
        for tlvIdx in range(numTLVs):
            
            # word array to convert 4 bytes to a 32 bit number
            word = [1, 2**8, 2**16, 2**24]
            
            # Check the header of the TLV message
            tlv_type = np.matmul(byteBuffer[idX:idX+4],word) # structure tag
            idX += 4
            tlv_length = np.matmul(byteBuffer[idX:idX+4],word) # length of structure
            idX += 4
            # Read the data depending on the TLV message
            range_depth = configParameters["numRangeBins"] * configParameters["rangeResolutionMeters"]
            range_width = range_depth / 2, 400
            if tlv_type == MMWDEMO_UART_MSG_DETECTED_POINTS:
                # word array to convert 2 bytes to a 16 bit number
                word = [1, 2**8]
                tlv_numObj = np.matmul(byteBuffer[idX:idX+2],word) 
                idX += 2 # ------------???
                tlv_xyzQFormat = 2**np.matmul(byteBuffer[idX:idX+2],word) # shouldnt it be 4 bytes after, not 2?
                idX += 2
                
                # Initialize the arrays
                rangeIdx = np.zeros(tlv_numObj,dtype = 'int16') # type int16 because they are 2 bytes values
                dopplerIdx = np.zeros(tlv_numObj,dtype = 'int16')
                peakVal = np.zeros(tlv_numObj,dtype = 'int16')
                x = np.zeros(tlv_numObj,dtype = 'int16')
                y = np.zeros(tlv_numObj,dtype = 'int16')
                z = np.zeros(tlv_numObj,dtype = 'int16')
                
                for objectNum in range(tlv_numObj):
                    
                    # Read the data for each object
                    rangeIdx[objectNum] =  np.matmul(byteBuffer[idX:idX+2],word)
                    idX += 2
                    dopplerIdx[objectNum] = np.matmul(byteBuffer[idX:idX+2],word)
                    idX += 2
                    peakVal[objectNum] = np.matmul(byteBuffer[idX:idX+2],word)
                    idX += 2
                    x[objectNum] = np.matmul(byteBuffer[idX:idX+2],word)
                    idX += 2
                    y[objectNum] = np.matmul(byteBuffer[idX:idX+2],word)
                    idX += 2
                    z[objectNum] = np.matmul(byteBuffer[idX:idX+2],word)
                    idX += 2
                    
                # Make the necessary corrections and calculate the rest of the data
                rangeVal = rangeIdx * configParameters["rangeIdxToMeters"]
                dopplerIdx[dopplerIdx > (configParameters["numDopplerBins"]/2 - 1)] = dopplerIdx[dopplerIdx > (configParameters["numDopplerBins"]/2 - 1)] - 65535
                dopplerVal = dopplerIdx * configParameters["dopplerResolutionMps"]
                x = x / tlv_xyzQFormat
                y = y / tlv_xyzQFormat
                z = z / tlv_xyzQFormat
                
                # Store the data in the detObj dictionary
                detObj = {"zi": [], "numObj": tlv_numObj, "rangeIdx": rangeIdx, "range": rangeVal, "dopplerIdx": dopplerIdx, \
                          "doppler": dopplerVal, "peakVal": peakVal, "x": x, "y": y, "z": z,"rangeArray":[], "rp":[], "np":[]}
                dataOK = 1  
            if tlv_type == MMWDEMO_UART_MSG_RANGE_PROFILE:
                numrp = 2 * configParameters["numRangeBins"]
                rp = byteBuffer[idX : idX + numrp]

                rp = list(map(add, rp[0:numrp:2], list(map(lambda x: 256 * x, rp[1:numrp:2]))))
                rp_x = (
                    np.array(range(configParameters["numRangeBins"]))
                    * configParameters["rangeIdxToMeters"]
                )
                idX += numrp
                detObj["rp"] = rp
                dataOK = 1      
            if tlv_type == MMWDEMO_UART_MSG_NOISE_PROFILE:
                numrp = 2 * configParameters["numRangeBins"]
                rp = byteBuffer[idX : idX + numrp]

                rp = list(map(add, rp[0:numrp:2], list(map(lambda x: 256 * x, rp[1:numrp:2]))))
                rp_x = (
                    np.array(range(configParameters["numRangeBins"]))
                    * configParameters["rangeIdxToMeters"]
                )
                idX += numrp
                detObj["np"] = rp
                dataOK = 1      
            if tlv_type == MMWDEMO_UART_MSG_RANGE_DOPPLER:
                numBytes = int(2*configParameters["numRangeBins"]*configParameters["numDopplerBins"])
                payload = byteBuffer[idX:idX + numBytes]
                idX += numBytes
                rangeDoppler = payload.view(dtype=np.int16)
            
                # Convert the range doppler array to a matrix
                rangeDoppler = np.reshape(rangeDoppler, (int(configParameters["numDopplerBins"]), configParameters["numRangeBins"]),'F') #Fortran-like reshape
                rangeDoppler = np.append(rangeDoppler[int(len(rangeDoppler)/2):], rangeDoppler[:int(len(rangeDoppler)/2)], axis=0)

                # Generate the range and doppler arrays for the plot
                rangeArray = np.array(range(configParameters["numRangeBins"]))*configParameters["rangeIdxToMeters"]
                dopplerArray = np.multiply(np.arange(-configParameters["numDopplerBins"]/2 , configParameters["numDopplerBins"]/2), configParameters["dopplerResolutionMps"])
                detObj = {"zi": [], "np": [], "rp":[], "rangeArray": rangeArray, "dopplerArray":dopplerArray, "rangeDoppler":rangeDoppler, "x":[]}
                dataOK = 1  
            if tlv_type == MMWDEMO_UART_MSG_AZIMUT_STATIC_HEAT_MAP:
                numTxAnt = 2
                numRxAnt = 4
                numBytes = numRxAnt * numTxAnt * configParameters["numRangeBins"] * 4
                q = byteBuffer[idX : idX + numBytes]
                idX += numBytes
                q = q[0::2] + q[1::2] * 2**8
                q[q > 32767] -= 65536
                q = q[0::2] + 1j * q[1::2]
                q = q.reshape((2*4, configParameters['numRangeBins']))
                Q = np.fft.fft(q, NUM_ANGLE_BINS)  # column based NUM_ANGLE_BINS-point fft, padded with zeros
                QQ = np.fft.fftshift(np.abs(Q), axes=0)
                QQ = QQ.T
                QQ = QQ[:, 1:]
                QQ = np.fliplr(QQ)
                theta = np.arcsin(np.arange(-NUM_ANGLE_BINS/2+1, NUM_ANGLE_BINS/2) * (2/NUM_ANGLE_BINS)) * (180/np.pi)
                range_ = np.arange(configParameters['numRangeBins']) * configParameters['rangeIdxToMeters']

                detObj["qq"] = QQ
                detObj["theta"] = theta
                detObj["range"] = range_
                dataOK = 1
            if tlv_type == MMWDEMO_UART_MSG_STATS:
                print("Getting called")
                print("STarting idX: ", idX)
                word = [1, 2**8, 2**16, 2**24]
                interFrameProcessingTime = np.matmul(byteBuffer[idX : idX + 4], word)
                idX += 4
                transmitOutputTime = np.matmul(byteBuffer[idX : idX + 4], word)
                idX += 4
                interFrameProcessingMargin = np.matmul(byteBuffer[idX : idX + 4], word)
                idX += 4
                interChirpProcessingMargin = np.matmul(byteBuffer[idX : idX + 4], word)
                idX += 4
                activeFrameCPULoad = np.matmul(byteBuffer[idX : idX + 4], word)
                idX += 4

                interFrameCPULoad = np.matmul(byteBuffer[idX : idX + 4], word)
                idX += 4

                statisticsObj = {
                    "interFrameProcessingTime": interFrameProcessingTime,
                    "transmitOutputTime": transmitOutputTime,
                    "interFrameProcessingMargin": interFrameProcessingMargin,
                    "interChirpProcessingMargin": interChirpProcessingMargin,
                    "activeFrameCPULoad": activeFrameCPULoad,
                    "interFrameCPULoad": interFrameCPULoad,
                }
                print(statisticsObj)
                print("Ending idX: ", idX)
        
  
        # Remove already processed data
        if idX > 0 and byteBufferLength > idX:
            shiftSize = totalPacketLen
               
            byteBuffer[:byteBufferLength - shiftSize] = byteBuffer[shiftSize:byteBufferLength]
            byteBuffer[byteBufferLength - shiftSize:] = np.zeros(len(byteBuffer[byteBufferLength - shiftSize:]),dtype = 'uint8')
            byteBufferLength = byteBufferLength - shiftSize
            
            # Check that there are no errors with the buffer length
            if byteBufferLength < 0:
                byteBufferLength = 0
                

    return dataOK, frameNumber, detObj

# ------------------------------------------------------------------

# Funtion to update the data and display in the plot
save_state_rp = []
save_state_np = []
saver = ""
def update():
    dataOk = 0
    global detObj, s, p, ptr, s2, p2, p3, s3
    x = []
    y = []
    peak_val = []
    range = []
    global save_state_rp, save_state_np,  noise_prof, saver
    # Read and parse the received data
    dataOk, frameNumber, detObj = readAndParseData14xx(Dataport, configParameters)
    if dataOk and "x" in detObj and len(detObj["x"]) > 0:
        x = -detObj["x"]
        y = detObj["y"]

        peak_val = detObj["peakVal"] # idk if those are dBm or what
        range = detObj["range"]

        v_doppler = detObj["doppler"]
        
        s.setData(x,y)
        s2.setData(range,peak_val)
        s3.setData(range,v_doppler)
        if ptr == 0:
            p.enableAutoRange('xy', False)  ## stop auto-scaling after the first data set is plotted
            p2.enableAutoRange('xy', False) 
        ptr += 1

    # if doppler heatmap received
    if dataOk and "rangeArray" in detObj and len(detObj["rangeArray"]) > 0:
        plt.title("Doppler Range Heatmap")
        plt.xlabel("Range (m)")
        plt.ylabel("Doppler (m/s)")
        cs = plt.contourf(detObj["rangeArray"],detObj["dopplerArray"],detObj["rangeDoppler"])
        fig.colorbar(cs, shrink=0.9)
        plt.pause(1e-10)
        plt.clf()
        # Append data to output CSV
        with open('output.csv', 'a+') as file:
                # Flatten Range-Doppler array from (16, 256) to (1, 4096) for easier saving
                csv.writer(file, delimiter=',').writerows([detObj["rangeArray"],detObj["dopplerArray"],list(chain.from_iterable(detObj["rangeDoppler"]))])
    
    #if azimuth heatmap received
    if dataOk and "qq" in detObj and len(detObj["qq"]) > 0:
        plt.title("Azimuth Range Heatmap")
        plt.xlabel("Range (m)")
        plt.ylabel("Doppler (m/s)")
        ax = plt.subplot(1, 1, 1)  # rows, cols, idx
        range_depth = configParameters["numRangeBins"] * configParameters["rangeResolutionMeters"]
        range_width, grid_res = range_depth / 2, 400
        plt.imshow(detObj["qq"], extent=[detObj["theta"][0], detObj["theta"][-1], detObj["range"][0], detObj["range"][-1]], vmin=0, vmax=np.max(detObj["qq"]), aspect='auto')
        ax = plt.gca()
        ax.invert_yaxis()
        ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 6)))
        plt.pause(1e-10)
        plt.clf()

    # if range profile received
    if dataOk and "rp" in detObj and len(detObj["rp"]) > 0:
        if not (save_state_np==[]):
            plt.title("Range/Noise Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            plt.plot(save_state_np, color='blue')
            peaks, _ = find_peaks(detObj["rp"], height=9500)
            plt.plot(peaks, [detObj["rp"][i] for i in peaks], "x")
            plt.plot(detObj["rp"], color='red')
            plt.pause(1e-10)
            plt.clf()
            # Append data to output CSV
            with open('output.csv', 'a+') as file:
                csv.writer(file, delimiter=',').writerows([detObj["rp"], save_state_np])
            save_state_np = []
        # if also tracking noise, save range profile and mark it as saver
        if noise_prof == 1:
            saver = "rp"
            save_state_rp = detObj["rp"]
        # if tracking range alone:
        else:
            plt.title("Range Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            peaks, _ = find_peaks(detObj["rp"], height=9500)
            cs = plt.plot(detObj["rp"], color='blue')
            plt.plot(peaks, [detObj["rp"][i] for i in peaks], "x")
            plt.pause(1e-10)
            plt.clf()
            # Append data to output CSV
            with open('output.csv', 'a+') as file:
                csv.writer(file, delimiter=',').writerows([detObj["rp"]])
    if dataOk and "np" in detObj and len(detObj["np"]) > 0:
        # if also tracking range and it has been saved beforehand
        if not (save_state_rp==[]):
            plt.title("Range/Noise Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            plt.plot(detObj["np"], color='red')
            peaks, _ = find_peaks(detObj["rp"], height=9500)
            plt.plot(peaks, [detObj["rp"][i] for i in peaks], "x")
            plt.plot(save_state_rp, color='blue')
            plt.pause(1e-10)
            plt.clf()
            # Append data to output CSV
            with open('output.csv', 'a+') as file:
                csv.writer(file, delimiter=',').writerows([save_state_rp, detObj["np"]])
            save_state_rp==[]
        # if also tracking range, save noise profile and mark it as saver
        if range_prof == 1:
            saver = "np"
            save_state_np = detObj["np"]
        # if tracking noise alone:
        else:
            plt.title("Noise Profile")
            plt.xlabel("Range (m)")
            plt.ylabel("Relative Power (dB)")
            ax = plt.gca()
            ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 100)))
            ax.xaxis.set_major_formatter(ticker.FuncFormatter(lambda x, pos: ('%g') % (x / 25)))
            cs = plt.plot(detObj["np"], color='red')
            plt.pause(1e-10)
            plt.clf()
            # Append data to output CSV
            with open('output.csv', 'a+') as file:
                csv.writer(file, delimiter=',').writerows([detObj["np"]])
    return dataOk


# -------------------------    MAIN   -----------------------------------------  


# Configurate the serial port
CLIport, Dataport = serialConfig(configFileName)

# Get the configuration parameters from the configuration file
configParameters = parseConfigFile(configFileName)

# Save GUI display type at top of CSV, create said file
if range_prof == 1:
    if noise_prof == 1:
        with open('output.csv', 'w') as file:
                file.write("range_noise")
                file.write('\n')
    else:
         with open('output.csv', 'w') as file:
                file.write("range")
                file.write('\n')
elif noise_prof == 1:
     with open('output.csv', 'w') as file:
                file.write("noise")
                file.write('\n')
if doppler_heat == 1:
     with open('output.csv', 'w') as file:
                file.write("doppler")
                file.write('\n')
ptr = 0

# START QtAPP for the plot
if det_obj:
    app = pg.mkQApp("Scattering Plot")

    win = pg.GraphicsLayoutWidget(show=True, title="Radar")
    win.resize(1000,600)
    win.setWindowTitle('pyqtgraph example: Plotting')
    pg.setConfigOptions(antialias=True)

    # Set the plot ----------------------------
    pg.setConfigOption('background','w')

    p = win.addPlot(title="Detected Objects: XY")
    p.setXRange(-1.5,1.5)
    p.setYRange(0,2)
    p.setLabel('left',text = 'Y position (m)')
    p.setLabel('bottom', text= 'X position (m)')
    p.showGrid(x=True, y=True, alpha=True)
    s = p.plot([],[],pen=None,symbol='x')

    win.nextRow()

    p2 = win.addPlot(title="Detected objects: Range")
    p2.setXRange(0,1.5)
    p2.setYRange(0,50)
    p2.setLabel('left',text = 'Peak value')
    p2.setLabel('bottom', text= 'Range (m)')
    p2.showGrid(x=True, y=True, alpha=True)
    s2 = p2.plot([],[],pen=None,symbol='x')

    win.nextRow()

    p3 = win.addPlot(title="Detected objects: Doppler velocity")
    p3.setXRange(0,1.5)
    p3.setYRange(0,2)
    p3.setLabel('left',text = 'Doppler velocity (m/s)')
    p3.setLabel('bottom', text= 'Range (m)')
    p3.showGrid(x=True, y=True, alpha=True)
    s3 = p3.plot([],[],pen=None,symbol='x')

    timer = QtCore.QTimer()
    timer.timeout.connect(update)
    timer.start(10)

    pg.exec()

# Main loop 
detObj = {}  
frameData = {}    
currentIndex = 0

fig = plt.figure()

while True:
    try:
        # Update the data and check if the data is okay
        dataOk = update()

        if dataOk:
            # Store the current frame into frameData
            frameData[currentIndex] = detObj
            currentIndex += 1

        # time.sleep(0.033) # Sampling frequency of 30 Hz
        
    # Stop the program and close everything if Ctrl + c is pressed
    except KeyboardInterrupt:
        CLIport.write(('sensorStop\n').encode())
        CLIport.close()
        Dataport.close()
        win.close()
        break
