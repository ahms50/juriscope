import socket
import numpy as np
import time
import csv
import pickle
import os
import re
import struct
from tifffile import TiffWriter
from datetime import datetime, timedelta, timezone
from ome_types.model import OME, Image, Pixels, Channel, Plane
import sys




################################## CHANGE EVERYTHING IN THIS SECTION ######################

load_previous_positions=1 # 1 to load previous positions you have already selected
                          # 0 to choose new positions for samples

main_folder='/home/ahm50/data/4_12_25/' # Give the file directory you are saving to

objective=40 # Set either 20 or 40 for which objective is being used

samples=1 #How many wells/capillaries do you have

sample_names=['Melting_in_GUV']

z_stack=1# Say if a z stack is being used; 1=yes, 0= not
timelapse=1 # say if you want time lapse on; 1=yes, 0= not

auto_focus=1 # Say if you want autofocus ON or OFF


z_c_order='zc'  # zc -  In single z-stack, does all channels then moves to next z_stack -- should be faster
                # cz -  In single c, does all z_stack and moves to next c --- should be more precise

##################################### time range

interval=np.array([5]) # interval time between imaging cycles in minutes


interval_time=np.array([1000]) #  total time for each interval, same order as interval

#interval_time=interval_time*60 # This is for interval_time in hrs.
                               # If using in minutes, comment out.

## Further intervals can be added

################################### z step range
z_range=5 # +- this number, e.g. if it is 20, it is -20Âµm to +20 Âµm
z_step=0.5 # step size between z range, i.e. if z_range is 20 and z_step is 1, you would have 41 slices

## Temperature and heating time

enclosure_heating=1 # Say if you want enclosure heating on; 1=yes 2=no

temperature_enclosure = 30
heating_time = "0:0:1" # h:m:s

#### Peltier temperature variables

peltier=1 # say if you want peltier on; 1=yes, 0=not

pelt_start_temp=25
pelt_end_temp=65
pelt_step=0.2

pelt_initial_heat_time=30 # how long the initial heating of the first temp lasts for before imaging begins in seconds

pelt_return=1 # Set: 0 -> Only goes from pelt_start_temp to pelt_end_temp in pelt_step
              #      1 -> Goes from pelt_start_temp to pelt_end_temp and goes back to pelt_start_temp in pelt_step


pelt_time_temp_threshold=[] # Cutoff temperature where the pelt_wait_time moves to the next wait time sequence
                                    # if only singular heating time interval used, let this be: []
                                    
pelt_wait_time=[5] # wait time for peltier in minutes. 
                              # SHOULD BE: length sshould be: len(pelt_time_temp_threshold) +1
                              # if only singular heating time interval, should only be 1 interval number and not list


################################## Illumination 


illumination=[0x40,0x20] # change depending on which channel is being used
illum_expose=[500000, 500000] # change depending on exposure time for each channel
laser_pwr=[1,1] # change depending on the laser power wanted, between 0 and 1


##Specify illumination setting
#first: define the illumination channel that you want to use:
#ch postion 1: BF 475 nm: 0 bit value , Hex value: 0x01
#ch postion 2: BF 528 nm: 2 bit value, Hex value: 0x02
#ch postion 3: BF 625 nm: 4 bit value, Hex value: 0x04
#ch postion 4: BF 655 nm: 8 bit value, Hex value: 0x08
#ch postion 5: Fluo 385 nm: 16 bit value, Hex value: 0x10
#ch postion 6: Fluo 475 nm: 32 bit value, Hex value: 0x20  --- error in microscope so turns on both 475 and 528 on the 40x
#ch postion 7: Fluo 528 nm: 64 bit value, Hex value: 0x40
#ch postion 8: Fluo 625 nm: 128 bit value, Hex value: 0x80


if load_previous_positions==1:
    #### Use this to load all positions data after an error. Make sure to comment out the bits above
    with open(f"{main_folder}/position_data.pkl", "rb") as f:
        number, positions = pickle.load(f)
        print("Loaded data:", number, positions)

###################################### END OF USER SET PARAMETERS #################################################



#################################### Initialisation

os.makedirs(main_folder, exist_ok=True)  # creates dircetory if missing


movie_output_path=f'{main_folder}/movies'

os.makedirs(movie_output_path, exist_ok=True)  # creates dircetory if missing

#assign cama name

if objective==20:
    camera_name = "Genicam FLIR Blackfly S BFS-U3-70S7M 0159F410" # 20x
elif objective==40:
    camera_name = "Genicam FLIR Blackfly S BFS-U3-70S7M 0159F411" #40x


########### Set peltier array

if pelt_return==0:
    pelt_temp_list=list(range(pelt_start_temp, pelt_end_temp - 1, pelt_step))
elif pelt_return==1:
    pelt_up=list(np.round(np.arange(pelt_start_temp, pelt_end_temp - 1, pelt_step),2))
    pelt_down=list(np.round(np.arange( pelt_end_temp, pelt_start_temp - 1, -pelt_step),2))
    pelt_temp_list = pelt_up + pelt_down


assert len(pelt_wait_time) == len(pelt_time_temp_threshold) + 1 # check if pelt wait time matches threshold

pelt_wait_schedule = []
for temp in pelt_temp_list:
    for i, threshold in enumerate(pelt_time_temp_threshold):
        if temp > threshold:
            pelt_wait_schedule.append(pelt_wait_time[i])
            break
    else:
        # If no threshold matched (i.e., temp <= all thresholds)
        pelt_wait_schedule.append(pelt_wait_time[-1])
        
        
########### Name mapping for illumination

# Channel mapping: hex value -> wavelength string
illumination_map = {
    0x01: "475 nm (BF)",
    0x02: "528 nm (BF)",
    0x04: "625 nm (BF)",
    0x08: "655 nm (BF)",
    0x10: "385 nm (Fluo)",
    0x20: "475 nm (Fluo)",
    0x40: "528 nm (Fluo)",
    0x80: "625 nm (Fluo)",
}

# Channel mapping laser power
illumination_val = {
    0x01: "0",
    0x02: "1",
    0x04: "2",
    0x08: "3",
    0x10: "4",
    0x20: "5",
    0x40: "6",
    0x80: "7",
}


# Translate hex values into wavelengths
illumination_names = [illumination_map[val] for val in illumination]





######################## Variables to edit

movie_folder_path=movie_output_path


output_path=f'{main_folder}/tiffs2'

os.makedirs(output_path, exist_ok=True)  # creates dircetory if missing

## Used this: illumination=[0x01,0x60,0x40,0x80] # change depending on which channel is being used

compression_type = "raw"  # or "tiff_lzw" for lossless compression or "raw" for no compresison

user_channel_colors = {
    #0: "#FFFFFFFF", # white for greyscale
    0: "#00FFFFFF", 
    1: "#FFFF00FF"  
    #3: "#FF0000FF"    
}


################################## variables and functions for movie2tiff conversions #################################


CAMERA_MOVIE_MAGIC = 0x496D6554  # 'TemI' little-endian
CAMERA_HEADER_LEN = 56 # character length of camera name

# defines the endians to check later
G_BIG_ENDIAN = 4321
G_LITTLE_ENDIAN = 1234





_STRUCT = struct.Struct("<7I2Q3I")  # Define Binary layout:
                                    # "<" : Little-endian byte order
                                    # 7I  : 7 unsigned 32-bit integers - all bits from FrameHeader magic -> endianness
                                    # 2Q  : 2 unsigned 64-bit integers - time_Sec and time_nsec
                                    # 3I  : 3 unsigned 32 bit integers - width, height, stride


def hdr_from_bytes(buf:bytes):
    if len(buf)<CAMERA_HEADER_LEN:
        raise ValueError("Incomplete camera_save_struct header")
    fields=_STRUCT.unpack_from(buf)


    # assign each of the unpacked binary information into their respective datasets
    hdr = {
        "magic": fields[0],
        "version": fields[1],
        "type": fields[2],
        "pixelformat": fields[3],
        "length_header": fields[4],
        "length_data": fields[5],
        "endianness": fields[6],
        "time_sec": fields[7],
        "time_nsec": fields[8],
        "width": fields[9],
        "height": fields[10],
        "stride": fields[11],
    }




    if hdr["magic"] != CAMERA_MOVIE_MAGIC: # makes sure that the buf starts with TemI
            raise ValueError("Bad TemI magic")
    

    return hdr



# Generator function to find the hdr: directory of fields and extra: data length of image frame

def iterate_frames(data):
    off=0
    total = len(data) # finds the total number of bytes of movie 

    magic=CAMERA_MOVIE_MAGIC.to_bytes(4,'little') # what TemI is in bytes

    while off<total: # looks through list upto the number of bytes in the movie
        idx=data.find(magic,off) # finds the index of the data which begins TemI starting from the off index

        if idx==-1:
            break

        extra_start=idx + CAMERA_HEADER_LEN # data index after the TemI file and camera name length are taken into account
        
        hdr=hdr_from_bytes(data[idx : extra_start]) # creates a dictionary from the function with relavent parts of the image
        
        
        extra_end= idx + hdr["length_header"] # where the bytes end for this frame
        extra=data[extra_start:extra_end] # looks at the movie file from indexes from extra_start to extra_end
        frame_start=extra_end # frame starts after all the heathers are finished
        frame_end= frame_start + hdr["length_data"] # all the data fromt the movie file for that specific frame
        yield hdr, extra, memoryview(data)[frame_start:frame_end] # memory view can access the bytes of data without copying them
        off=frame_end # Start of the next iteration



################################## Functions for ome.tiff conversion ####################



# -----------------------
# OMEWriter class 
# -----------------------
class OMEWriter:
    def __init__(self, output_path, filename, shape, timestamp_list, user_channel_colors):
        T, C, Z, Y, X = shape

        PREDEFINED_COLORS = {
            "grey":    '#808080FF',
            "red":     '#FF0000FF',
            "green":   '#00FF00FF',
            "blue":    '#0000FFFF',
            "yellow":  '#FFFF00FF',
            "cyan":    '#00FFFFFF',
            "magenta": '#FF00FFFF',
            "white":   '#FFFFFFFF',
            "black":   '#000000FF',
        }

        def hex_to_argb(hex_str):
            hex_str = hex_str.lstrip("#")
            if len(hex_str) == 6:
                argb = 0xFF000000 | int(hex_str, 16)
            elif len(hex_str) == 8:
                argb = int(hex_str, 16)
            else:
                raise ValueError("Invalid color")
            return argb if argb < (1 << 31) else argb - (1 << 32)

        # Build OME Channels
        channels = []
        for i in range(C):
            col = user_channel_colors.get(i, "white")
            if not col.startswith("#"):
                col = PREDEFINED_COLORS[col]
            col = hex_to_argb(col)
            channels.append(Channel(id=f"Channel:{i}", name=f"C{i}", samples_per_pixel=1, color=col))

        # Build Planes
        planes = []
        for t in range(T):
            for c in range(C):
                for z in range(Z):
                    planes.append(
                        Plane(
                            the_t=t,
                            the_c=c,
                            the_z=z,
                            delta_t=float(timestamp_list[t])
                        )
                    )

        pixels = Pixels(
            dimension_order="XYZCT",
            type="uint16",
            size_x=X, size_y=Y,
            size_z=Z,
            size_c=C,
            size_t=T,
            channels=channels,
            planes=planes
        )

        ome = OME(images=[Image(id="Image:0", name=filename, pixels=pixels)])
        self.ome_xml = ome.to_xml()

        self.tif = TiffWriter(f"{output_path}/{filename}.ome.tif", bigtiff=True)
        self.is_first = True

    def write_plane(self, arr):
        """Write one (Y,X) plane."""
        self.tif.write(
            arr,
            compression="lzw",
            photometric="minisblack",
            description=self.ome_xml if self.is_first else None
        )
        self.is_first = False

    def close(self):
        self.tif.close()


# List all files in the folder
files = [f for f in os.listdir(movie_folder_path) if os.path.isfile(os.path.join(movie_folder_path, f)) and f.lower().endswith(".movie")]



############## Obtain the different samples from this.
# Dictionary to store files grouped by identifier
grouped_files = {}

if enclosure_heating==1:

    for file in files: #read every single files
        if file.startswith("sample_") and "_enclosure" in file:
            # Extract the identifier between 'sample_' and '_enclosure'
            identifier = file.split("sample_")[1].split("_enclosure")[0]
            
            if identifier not in grouped_files:
                grouped_files[identifier] = []
            
            grouped_files[identifier].append(file)

        # print(grouped_files)
else:
    if peltier==1:
        for file in files: #read every single files
            if file.startswith("sample_") and "_peltier" in file:
                # Extract the identifier between 'sample_' and '_enclosure'
                identifier = file.split("sample_")[1].split("_peltier")[0]
                
                if identifier not in grouped_files:
                    grouped_files[identifier] = []
                
                grouped_files[identifier].append(file)

            # print(grouped_files)
    
    
    
    else:
    
    
        for file in files: #read every single files
            if file.startswith("sample_") and "." in file:
                # Extract the identifier between 'sample_' and '_enclosure'
                identifier = file.split("sample_")[1].split(".")[0]
                
                if identifier not in grouped_files:
                    grouped_files[identifier] = []
                
                grouped_files[identifier].append(file)

            # print(grouped_files)
        


############ Sort each group in increasing timestamp output is grouped_files
for identifier, files_list in grouped_files.items():
    # Extract 'n' and sort
    files_list.sort(key=lambda f: int(re.search(r'_timestamp_(\d+)_', f).group(1)))


sample_names=list(grouped_files.keys()) # obtain sample names



########### open all frame mapping which contains the data:
                            #"frame": frame_seq_id,
                            #"sample": sample_names,
                            #"number": num,
                            #"z": zi,
                            #"illum_wavelength": illumination,
                            #"illum_exposure_time": illum_expose,
                            #"illum_pwr": laser_pwr


########### Have them in all_data, split between different samples
all_data = []  # This will hold multiple instances of data

for sample_idx in range(len(sample_names)):
     
     # Open the CSV file
    with open(f"{main_folder}/frame_mapping_{sample_idx}.csv", mode='r', newline='') as file:
        reader = csv.DictReader(file)
        data = [row for row in reader]  # Each row is a dictionary
        all_data.append(data)  # Store this instance of data


########## CAN CREATE SEQUENCE ARRAY USING FRAME NUMBERS AND CROSS REFERENCE WITH TIMESTAMPS. 
########## THEN ITERATE TIMESTAMPS AS PER WHAT IS GIVEN. - MAKE SURE TO CHANGE BETWEEN FILES. CAN BE DONE BY ITERATING TIMESTAMP WITH FILE ORDER, AS IT IS ORDERED

#################### iterate through frames 


for sample_idx in range(len(sample_names)):

    num_pos = sorted(set(d.get('number') for d in all_data[sample_idx] if 'number' in d)) # find number of positions per sample

    for num in num_pos:  # for number of samples in the movie files

        num_dicts = [d for d in all_data[sample_idx] if d.get('number')==num] # takes only a specific number from the sample
        c_per_num=sorted(set(d.get('illum_wavelength') for d in num_dicts if 'illum_wavelength' in d)) # find number of channels per sample

        Y, X = 2200, 3208

        total_t=len(grouped_files[next(iter(grouped_files))])
        total_c=len(set(d.get('illum_wavelength') for d in num_dicts if 'illum_wavelength' in d))
        total_z=len(set(d.get('z') for d in num_dicts if 'z' in d))


        timestamp_list = []

        for ti, file in enumerate(grouped_files[sample_names[sample_idx]]):
            # Extract timestamp from filename
            timestamp_list.append(int(re.search(r'_timestamp_(\d+)_', file).group(1)))


        # --- Create writer ---
        writer = OMEWriter(
            output_path=output_path,
            filename=f"sample_{sample_names[sample_idx]}_position_{num}",
            shape=(total_t, total_c, total_z, Y, X),
            timestamp_list=timestamp_list,
            user_channel_colors=user_channel_colors
        )
     
              

        for ti, file in enumerate(grouped_files[sample_names[sample_idx]]): # takes the ordered timelapse for each sample and reads them in order of timelapse T
            movie_name=f"{movie_folder_path}/{file}"
            with open(movie_name, "rb") as f:
                data = f.read() # reads the movie file as bytes

            timestamp_list.append(int(re.search(r'_timestamp_(\d+)_', file).group(1))) 
            
            frames=list(iterate_frames(data)) # creates a list of all frames from data
            n_frames=len(frames) # finds the number of frames

            if n_frames == 0:
                    raise ValueError("No TemI frames found") # incase no frames are in image

            #print(n_frames)
            #print(file)
            
    
        
            for ci,c in enumerate(c_per_num): # for number of channels in sample
                c_dicts=[d for d in num_dicts if d.get('illum_wavelength')==c] # takes only a specific channel from the sample
                
                z_per_c=sorted({float(d['z']) for d in c_dicts if 'z' in d})
                
                for zi,z in enumerate(z_per_c): # for number of z_stacks in channel
                    z_dicts=[d for d in c_dicts if d.get('z')==str(z)]
                    
                    frame=set(d.get('frame') for d in z_dicts if 'frame' in d)
                    frame = int(list(frame)[0])

                    if len(all_data[sample_idx]) + 1 == n_frames: # if statement as sometimes 1st frame is incorrect for use
                        frame_idx=frame+1
                    elif len(all_data[sample_idx]) == n_frames:
                        frame_idx=frame
                    else:
                        raise ValueError("Frames do not match")
                    
                    ############# Read and store frames in OME-Tiff
                    
                    hdr, extra, mv = frames[frame_idx - 1] # indexing for individual frame

                    h , w = hdr["height"] , hdr["width"] # finds the height and width of the frame
                    stride=hdr["stride"] # difference in bytes between start of 1 row and the next row of pixels

                    arr=np.empty((h, w), dtype=np.uint16) # output image is of the 16 bit type. This initialises the matrix
                    for r in range(h):
                        off = r * stride # number of bytes in each row, i.e. total number of bytes in that row
                        row=np.frombuffer(mv[off : off + w * 2], dtype=np.uint16, count=w) # creates array from buffer of binary data using memory view between off: off + w and reads w number of elements 
                        if hdr["endianness"] == G_BIG_ENDIAN: # check endianess and make sure it is little endian, i.e. byte order is 12 34 and not 43 21
                            row = row.byteswap()
                        arr[r] = row
                    
                    #frame_time=hdr["time_sec"] + hdr["time_nsec"]/1e9
                    #iso_time = datetime.fromtimestamp(frame_time, tz=timezone.utc).isoformat().replace("+00:00", "Z")

                    #delta_t = timestamp_list[ti]




                    writer.write_plane(arr)
        writer.close()
                    

print('Conversion to tiff complete')

