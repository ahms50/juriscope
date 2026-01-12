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

main_folder='/home/ahm50/data/8_1_26/' # Give the file directory you are saving to

objective=40 # Set either 20 or 40 for which objective is being used

samples=2 #How many wells/capillaries do you have

sample_names=['Melting_Bulk_no_anch','Melting_Bulk_2uM_anch']

z_stack=1# Say if a z stack is being used; 1=yes, 0= not
timelapse=1 # say if you want time lapse on; 1=yes, 0= not

auto_focus=0 # Say if you want autofocus ON or OFF


z_c_order='zc'  # zc -  In single z-stack, does all channels then moves to next z_stack -- should be faster
                # cz -  In single c, does all z_stack and moves to next c --- should be more precise

##################################### time range

interval=np.array([5]) # interval time between imaging cycles in minutes


interval_time=np.array([1500]) #  total time for each interval, same order as interval

#interval_time=interval_time*60 # This is for interval_time in hrs.
                               # If using in minutes, comment out.

## Further intervals can be added

################################### z step range
z_range=10 # +- this number, e.g. if it is 20, it is -20Âµm to +20 Âµm
z_step=2 # step size between z range, i.e. if z_range is 20 and z_step is 1, you would have 41 slices

## Temperature and heating time

enclosure_heating=1 # Say if you want enclosure heating on; 1=yes 2=no

temperature_enclosure = 30
heating_time = "0:0:1" # h:m:s

#### Peltier temperature variables

peltier=1 # say if you want peltier on; 1=yes, 0=not

pelt_start_temp=25
pelt_end_temp=55
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


illumination=[0x20] # change depending on which channel is being used
illum_expose=[10000] # change depending on exposure time for each channel
laser_pwr=[1] # change depending on the laser power wanted, between 0 and 1


# remember to not have a comma for the last one
user_channel_colors = {
    #0: "#FFFFFFFF", # white for greyscale
    0: "#00FFFFFF"
    #1: "#FFFF00FF"
    
    #3: "#FF0000FF"    
}


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

############# initialise timestamp list
timestamp=[]


##### Open and configure a socket for allow the comunication between the py script and the temika software 
#the Host and Port are the same for both microscopes
HOST = "127.0.0.1"
PORT = 60000
sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM) # AF_INET: set address family as IPv4, SOCK_STREAM sockeet type TCP
sock.connect((HOST, PORT)) # connect to the host and the port given
sock.settimeout(120.) #timeout for connection is set at 1 second
reply = sock.recv(1024)

new_pos_save_vec=[]



####################################### Funtions ####################################




def script_send_command(command, reply_wait, terminator=b'\n\r', max_wait=120):
    """
    Sends a command to the socket and optionally waits for a reply ending in `terminator`.

    Args:
        command (str or bytes): The command to send.
        reply_wait (bool): Whether to wait for a reply.
        terminator (bytes): Sequence marking the end of a message.
        max_wait (int): Maximum time to wait for reply in seconds.

    Returns:
        bytes or None: Reply without terminator, or None if no reply expected.
    """

    if isinstance(command, str):
        command = command.encode('utf-8')

    sock.sendall(command)

    if reply_wait:
        reply = b""
        sock_start_time = time.time()

        while True:
            if time.time() - sock_start_time > max_wait:
                raise TimeoutError(f"No full reply received within {max_wait} seconds.")

            try:
                chunk = sock.recv(1024)
                if not chunk:
                    raise ConnectionError("Socket connection closed by the remote host.")

                reply += chunk

                if reply.endswith(terminator):
                    return reply[:-len(terminator)]

            except socket.timeout:
                continue  # just keep waiting

    else:
        time.sleep(0.1)


def wait_for_stepper_status(axis: str) -> str:
    """
    Waits for the stepper axis and status replies for the given axis,
    then returns the reply string once it's valid and contains a truthy float.
    
    Args:
        axis (str): The axis to query, e.g. 'x', 'y', or 'z'.
    
    Returns:
        str: The decoded status reply.
    """

    axis_command = f'<stepper axis="{axis}"></stepper>'.encode('utf-8')
    expected_axis_reply = f'MICROSCOPEONE STEPPER axis="{axis}"'

    # Step 1: Wait for correct axis reply
    while True:
        reply = script_send_command(axis_command, True)

        decoded_reply = reply.decode('utf-8').strip()



        if decoded_reply == expected_axis_reply:
            break


    # Send stepper axis command without expecting immediate reply
    script_send_command(f'<stepper axis="{axis}">'.encode('utf-8'), False)

    # Step 2: Wait for valid status reply
    while True:
        reply = script_send_command(b'<status></status>', True)


        decoded_reply = reply.decode('utf-8').strip()
        words = decoded_reply.split()

        if len(words) >= 4 and words[0] == 'MICROSCOPEONE' and words[1] == 'STEPPER' and words[2] == 'status':
            try:
                value = float(words[3])
                if value:
                    break
            except ValueError:
                # Handle the case where words[3] isn't a valid float
                pass



    # Send closing stepper tag
    script_send_command(b'</stepper>', False)

    # Return the full decoded reply for further processing
    return decoded_reply


# Function to get gets the current x, y, z, and a stepper motor positions from the microscope
def script_get_position():
    script_send_command(b'<microscopeone>', False) # sends command to open <microscopeone>


    reply = wait_for_stepper_status('x')

    
    x = float(reply.split()[3]) # receives the data anbd splits it where there is a whitespace and and takes the 4th word. In this case it is the x value
    
    
    # process repeated for y,z and apprature
    reply = wait_for_stepper_status('y')

    y = float(reply.split()[3])

    reply = wait_for_stepper_status('z')

    z = float(reply.split()[3])

    reply = wait_for_stepper_status('a')

    a = float(reply.split()[3])

    script_send_command(b'</microscopeone>', False)

    return [x, y, z, a]

# Function to tell the microscope where to move without pfs
def script_print_move_no_pfs(vec):
    script_send_command("\t<microscopeone>\n", False)
    script_send_command("\t\t<afocus>\n",False)
    script_send_command("\t\t\t<enable>OFF</enable>\n",True)
    script_send_command("\t\t</afocus>\n", False)

    axes = ['x', 'y', 'z']
    speeds = [10000, 10000, 100]

    for axis, value, speed in zip(axes, vec, speeds): #The zip() function takes multiple iterables (like lists) as arguments and aggregates elements from each of them.
        script_send_command(f'\t\t<stepper axis="{axis}">\n',True)
        script_send_command(f"\t\t\t<move_absolute>{value} {speed}</move_absolute>\n",False)
        #time.sleep(1)
        script_send_command("\t\t</stepper>\n", False)

    for axis in axes:
        script_send_command(f'\t\t<stepper axis="{axis}">\n', True)
        script_send_command("\t\t\t<wait_moving_end></wait_moving_end>\n", True)
        script_send_command("\t\t</stepper>\n\n", False)

    script_send_command("\t</microscopeone>\n", False)

# Function to tell the microscope where to move with pfs on
def script_print_move(vec):
    script_send_command("\t<microscopeone>\n", False)

    script_send_command("\t\t<afocus>\n", False)
    script_send_command("\t\t\t<enable>OFF</enable>\n", True)
    script_send_command("\t\t</afocus>\n", False)

    axes = ['x', 'y', 'z', 'a']
    speeds = [10000, 10000, 100, 10000]

    for axis, value, speed in zip(axes, vec, speeds):
        script_send_command(f'\t\t<stepper axis="{axis}">\n', True)
        script_send_command(f"\t\t\t<move_absolute>{value} {speed}</move_absolute>\n", False)
        #time.sleep(3) # change to a bit higher for cappillaries
        script_send_command("\t\t</stepper>\n", False)

    for axis in axes:
        script_send_command(f'\t\t<stepper axis="{axis}">\n', True)
        script_send_command("\t\t\t<wait_moving_end></wait_moving_end>\n", True)
        script_send_command("\t\t</stepper>\n\n", False)

    if auto_focus==1:
        script_send_command("\t\t<afocus>\n", False)
        script_send_command("\t\t\t<enable>ON</enable>\n", True)
        script_send_command("\t\t\t<wait_lock>0.2 10.3</wait_lock>\n", True) 
        script_send_command("\t\t</afocus>\n", False)
    elif auto_focus==0:
        script_send_command("\t\t<afocus>\n", False)
        script_send_command("\t\t\t<enable>OFF</enable>\n", True)
        script_send_command("\t\t</afocus>\n", False)

    script_send_command("\t</microscopeone>\n", False)


def peltier_temp(pelt_temp,pelt_wait):
    script_send_command("\t<microscopeone>\n", False)
    script_send_command('\t\t<pwr number="0">\n', False)
    script_send_command("\t\t\t<feedback>\n", False)
    script_send_command(f"\t\t\t\t<set>{pelt_temp}</set>\n", True)
    script_send_command("\t\t\t\t<enable>ON</enable>\n", True)
    script_send_command("\t\t\t</feedback>\n", False)
    script_send_command("\t\t</pwr>\n", False)
    script_send_command("\t</microscopeone>\n", False)
    
    time.sleep(pelt_wait)



# Function to tell the microscope to go in resting mode, pull fully down the objective and fully up the condeser, no pfs
def script_print_move_z_c_reset():
    script_send_command("\t<microscopeone>\n", False)
    
    #disable the pfs
    script_send_command("\t\t<afocus>\n", False)
    script_send_command("\t\t\t<enable>OFF</enable>\n", True)
    script_send_command("\t\t</afocus>\n", False)

    #move the objective (z axis) to -10000 and woth speed 100
    script_send_command(f'\t\t<stepper axis="z">\n', False)
    script_send_command(f"\t\t\t<move_absolute>-10000 100</move_absolute>\n", False)
    #time.sleep(100)
    script_send_command("\t\t</stepper>\n", False)

    script_send_command(f'\t\t<stepper axis="z">\n', True)
    script_send_command("\t\t\t<wait_moving_end></wait_moving_end>\n", True)
    script_send_command("\t\t</stepper>\n\n", False)
     
    #move the condenser (c axis) to +720000 and woth speed 30000
    script_send_command(f'\t\t<stepper axis="c">\n', False)
    script_send_command(f"\t\t\t<move_absolute>72000 30000</move_absolute>\n", False)
    #time.sleep(3)
    script_send_command("\t\t</stepper>\n", False)

    script_send_command(f'\t\t<stepper axis="c">\n', False)
    script_send_command("\t\t\t<wait_moving_end></wait_moving_end>\n", True)
    script_send_command("\t\t</stepper>\n\n", False)

    script_send_command("\t</microscopeone>\n", False)


# function for vectors to obtain z stack
def script_z_stack(vec,z_range,step):
    x, y, z, a = vec
    return [[x, y, z + offset, a] for offset in range(-z_range, z_range + 1, step)]

# Function for vectors to obtain z stack
def script_print_move_z_move(step):
    script_send_command("\t<microscopeone>\n",False)
 
    script_send_command("\t\t<afocus>\n",False)
    script_send_command("\t\t\t<enable>OFF</enable>\n",True)
    script_send_command("\t\t</afocus>\n",False)
 
    script_send_command(f'\t\t<stepper axis="z">\n',True)
    script_send_command(f"\t\t\t<move_absolute>{step} 100</move_absolute>\n",True)
    
    script_send_command("\t\t</stepper>\n",False)
 
    script_send_command(f'\t\t<stepper axis="z">\n',True)
    script_send_command("\t\t\t<wait_moving_end></wait_moving_end>\n",True)
    script_send_command("\t\t</stepper>\n\n",True)
 
    script_send_command("\t</microscopeone>\n",False)

### Main function which enables does the recording 
def script_print_sample(filename, positions, number):

    global new_pos_save_vec
    #Writes XML commands for imaging all positions in a capillary/well, including metadata, camera/microscope control, and image acquisition.
    #saves the metadata 
    #script_send_command("<!-- Image all positions in a capillary/well with z -->\n", False)
    #script_send_command("\t<!-- Start recording. -->\n", False)
    script_send_command("\t<save>\n", False)
    script_send_command(f"\t\t<header>\n{filename}\nNumber of positions {number}.\n\t\t</header>\n", True)
    script_send_command(f"\t\t<basename>{movie_output_path}/{filename}</basename>\n", True) # This is the file_location
    script_send_command("\t\t<append>DATE</append>\n", True)
    script_send_command("\t</save>\n", False)

    ########### acquisition

    ## Turn on camera with trigger mode
    script_send_command(f'\t<camera name="{camera_name}">\n', False)
    script_send_command(f'\t<genicam><enumeration feature="TriggerMode">On</enumeration></genicam>\n', False)
    script_send_command("\t\t<record>ON</record>\n", True)
    script_send_command("\t</camera>\n", False)

    ## Turn on metadata for environment
    if peltier == 1 and timelapse == 1:
        script_send_command("\t<microscopeone>\n", False)
        script_send_command("\t\t<record>ON</record>\n", True)
        script_send_command("\t</microscopeone>\n\n", False)

    ## Move and acquire image
    #script_send_command("\t<!-- Move to a position and acquire images -->\n", False)

    for i in range(number):

        if z_stack == 1:
            if z_c_order=='zc':

                pos_vec=positions[i] # select vector in position i
                script_print_move(pos_vec) # moves to the given vector

                new_vec=script_get_position()

                new_vec_time=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
                new_pos_time=[new_vec_time]+new_vec
                new_pos_save_vec.append(new_pos_time)

                
                for step in np.arange(-z_range,z_range +z_step,z_step):
                    script_print_move_z_move(new_vec[2]+step)

                    for illum, exp, laser in zip(
                        illumination,
                        illum_expose,
                        laser_pwr
                    ):
                        
                        #records the data given the illumination and exposure time
                        script_send_command("\t<microscopeone>\n", False)
                        script_send_command("\t\t<illumination>\n", False)
                        script_send_command(f"\t\t\t<enable>{illum}</enable>\n", True)
                        script_send_command(f'\t\t\t<value number="{illumination_val[illum]}">{laser}</value>\n', True) # double check number
                        script_send_command("\t\t</illumination>\n", False)
                        script_send_command("\t</microscopeone>\n", False)
                        script_send_command(f'\t<camera name="{camera_name}">\n', False)
                        script_send_command(f'\t<genicam><float feature="ExposureTime">{exp}</float></genicam>\n', True)
                        script_send_command("\t\t<send_trigger></send_trigger>\n", False)
                        script_send_command("\t</camera>\n\n", True)
                        # Turn off laser power for all channels
                        script_send_command("\t<microscopeone>\n", False)
                        script_send_command("\t\t<illumination>\n", False)
                        script_send_command(f"\t\t\t<enable>0</enable>\n", False)
                        script_send_command("\t\t</illumination>\n", False)
                        script_send_command("\t</microscopeone>\n", False)

                        if timelapse ==1: 

                            image_time=time.time()
                            timestamp_time=(image_time-overall_start_time) / 60

                            timestamp.append(timestamp_time)

            elif z_c_order=='cz':

                pos_vec=positions[i] # select vector in position i
                script_print_move(pos_vec) # moves to the given vector

                new_vec=script_get_position()

                new_vec_time=datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        
                new_pos_time=[new_vec_time]+new_vec
                new_pos_save_vec.append(new_pos_time)

                for illum, exp, laser in zip(
                        illumination,
                        illum_expose,
                        laser_pwr
                    ):

                    
                    
                    for step in np.arange(-z_range,z_range +z_step,z_step):
                        script_print_move_z_move(new_vec[2]+step)

                        #records the data given the illumination and exposure time
                        script_send_command("\t<microscopeone>\n", False)
                        script_send_command("\t\t<illumination>\n", False)
                        script_send_command(f"\t\t\t<enable>{illum}</enable>\n", True)
                        script_send_command(f'\t\t\t<value number="{illumination_val[illum]}">{laser}</value>\n', True) # double check number
                        script_send_command("\t\t</illumination>\n", False)
                        script_send_command("\t</microscopeone>\n", False)
                        script_send_command(f'\t<camera name="{camera_name}">\n', False)
                        script_send_command(f'\t<genicam><float feature="ExposureTime">{exp}</float></genicam>\n', True)
                        script_send_command("\t\t<send_trigger></send_trigger>\n", False)
                        script_send_command("\t</camera>\n\n", True)
                        # Turn off laser power for all channels
                        script_send_command("\t<microscopeone>\n", False)
                        script_send_command("\t\t<illumination>\n", False)
                        script_send_command(f"\t\t\t<enable>0</enable>\n", False)
                        script_send_command("\t\t</illumination>\n", False)
                        script_send_command("\t</microscopeone>\n", False)

                        if timelapse==1:

                            image_time=time.time()
                            timestamp_time=(image_time-overall_start_time) / 60

                            timestamp.append(timestamp_time)
                        
                        

        else:
                vec = positions[i]
                script_print_move(vec)

                for illum, exp, laser in zip(illumination, illum_expose, laser_pwr):
                    script_send_command("\t<microscopeone>\n", False)
                    script_send_command("\t\t<illumination>\n", False)
                    script_send_command(f"\t\t\t<enable>{illum}</enable>\n", False)
                    script_send_command(f'\t\t\t<value number="{illumination_val[illum]}">{laser}</value>\n', True)
                    script_send_command("\t\t</illumination>\n", False)
                    script_send_command("\t</microscopeone>\n", False)
                    script_send_command(f'\t<camera name="{camera_name}">\n', False)
                    script_send_command(f'\t<genicam><float feature="ExposureTime">{exp}</float></genicam>\n', True)
                    script_send_command("\t\t<send_trigger></send_trigger>\n", False)
                    script_send_command("\t</camera>\n\n", False)
                    # Turn off laser power for all channels
                    script_send_command("\t<microscopeone>\n", False)
                    script_send_command("\t\t<illumination>\n", False)
                    script_send_command(f"\t\t\t<enable>0</enable>\n", False)
                    script_send_command("\t\t</illumination>\n", False)
                    script_send_command("\t</microscopeone>\n", False)

                    if timelapse ==1:

                        image_time=time.time()
                        timestamp_time=(image_time-overall_start_time) / 60

                        timestamp.append(timestamp_time)
                        
    # Stop the recording for this range of values
    script_send_command("\t<!-- Stop recording. -->\n", False)
    script_send_command(f'\t<camera name="{camera_name}">\n', False)
    script_send_command(f'\t\t<genicam><enumeration feature="TriggerMode">Off</enumeration></genicam>\n', True) # turn off trigger mode
    script_send_command("\t\t<record>OFF</record>\n", True)
    script_send_command("\t</camera>\n", False)

    # metadata record off
    if peltier == 1 and timelapse == 1:
        script_send_command("\t<microscopeone>\n", False)
        script_send_command("\t\t<record>OFF</record>\n", True)
        script_send_command("\t</microscopeone>\n\n", False)
	



################################### get positions from temika ###################################

script_send_command(b"<temika>", False) #n.b. keep temika open at all times otherwise it cannot close


## If loading positions data, comment out from here

if load_previous_positions==0:
    print("p -- Save position\nu -- undo last capillary\nf -- Finish capillary\n")
    i = 0 # counter for how many capillaries/wells I am looking at
    number = [0] # number of positions in a given capillary/well
    positions = [] # stores all positions as a list, every 4 numbers are a position set
    while i < samples:
        key = input() # wait for the user to press a key

        # save the current position
        if key == 'p':
            vec = script_get_position() # obtains current position of the microscope
            number[i] += 1 # for a capillary/well, how many positions are being chosen
            positions.append(vec) # saves the positions
            print("capillary/well", i, vec) # prints what has been chosen

        #undo last saved position
        elif key == 'u':
            if number[i] > 0:
                number[i] -= 1 
                removed_vec=positions.pop() #removes last vector which was saved
                print("Last position removed from capillary/well",i,removed_vec)

        #finalise positions for the current well
        elif key == 'f':
            i += 1
            if i<samples:
                number.append(0)

            print("p -- Save position\nf -- Finish capillary\n")


    with open(f"{main_folder}/position_data.pkl", "wb") as f:
        pickle.dump((number, positions), f)
    print("Positions saved")

elif load_previous_positions==1:
    #### Use this to load all positions data after an error. Make sure to comment out the bits above
    with open(f"{main_folder}/position_data.pkl", "rb") as f:
        number, positions = pickle.load(f)
        print("Loaded data:", number, positions)



#print all the capillaries/wells, their positions and the vectors
i = 0
index = 0
while i < samples:
    print("capillary/well ", i)
    for j in range(number[i]):
        print("    position ", j, positions[j + index])
    index += number[i]
    i += 1




######################### check if everything is correct before proceeding

print(f'Check all settings')

print(f'Save folder: {main_folder}')

print(f'Sample names: {sample_names}')

print(f'Auto focus: {'ON' if auto_focus else 'OFF'}')

print(f'Channels and order: {', '.join(illumination_names)}')
print(f'Illumination Power: {laser_pwr}')
print(f'Exposure time: {illum_expose}')

print(f"z stack : {'Yes' if z_stack else 'No'}")

if z_stack==1:
    print(f'z range = {z_range}')
    print(f'z step = {z_step}')
    
    print(f'zc order: {z_c_order}')

print(f'Peltier = : {'Yes' if peltier else 'No'}"')

if peltier==1:
    print(f'Pelt start temp = {pelt_start_temp}')
    print(f'Pelt end temp = {pelt_end_temp}')
    print(f'Pelt temp step = {pelt_step}')
    
    print('Pelt temp array =', [f'{x:.2f}' for x in pelt_temp_list])
    print(f'Cut-off temperatures for waiting are: {pelt_time_temp_threshold}')
    print(f'Interval time per cut-off temperature range: {pelt_wait_time}')
    
    print(f'Heating time to reach 1st temperature: {pelt_initial_heat_time}s')

if enclosure_heating==1:
    print(f'Enclosure temp = {temperature_enclosure}')

print(f'timelapse = : {'Yes' if timelapse else 'No'}')

if timelapse==1:
    print(f'time interval (mins) = {interval}')
    print(f'Duration per interval (mins) = {interval_time}')
    print(f'Duration per interval (hrs) = {interval_time/60}')
    
    

total_time=sum(interval_time)

print(f'Total time (mins) = {total_time}')


print("y - continue, n - Terminate")

key = input() # wait for the user to press a key

if key == 'y':
    print("Continuing...")
    # your continuation code here
elif key == 'n':
    print("Terminating.")
    # optionally exit the script
    script_send_command(b"</temika>", False) #n.b. keep temika open at all times otherwise it cannot close
    exit()

################################## Start the script
#script_send_command(b"<temika>", False)
# Set the enclosure temperature and wait for stabilization
script_send_command(b"<!-- Set temperature_enclosure_high and wait heating_time_enclosure -->\n", False)

if enclosure_heating==1:

    script_send_command(b"<microscopeone>", False)

    script_send_command(b'\t\t<pwr number="4">\n', True) # 0: peltier, 1:PWR_0, 2:PWR_1 3: PWR_2, 4:encloser, 5:LED, 6:Fan
    script_send_command(b"\t\t\t<feedback>\n", False)
    script_send_command(f"\t\t\t\t<set>{temperature_enclosure}</set>\n", True)
    script_send_command(b"\t\t\t\t<enable>ON</enable>\n", True)
    script_send_command(b"\t\t\t</feedback>\n", False)
    script_send_command(b"\t\t</pwr>\n", False)
    script_send_command(b'<pwr number="5"><duty_cycle>0</duty_cycle></pwr>', True)
    script_send_command(b"\t</microscopeone>\n", False)



    h, m, s = map(int, heating_time.split(":"))
    heating_time_s = h * 3600 + m * 60 + s

    print(f"Enclosure heating started for {h:02}:{m:02}:{s:02}")

    # Countdown loop
    for remaining in range(heating_time_s, 0, -1):
        mins, secs = divmod(remaining, 60)
        hours, mins = divmod(mins, 60)
        print(f"\rTime remaining: {hours:02}:{mins:02}:{secs:02}", end="")
        time.sleep(1)
    print("\rTime remaining: 00:00:00")

    print("\nHeating complete. Imaging Start.")
    print('------------------------------------------------------------------')


if timelapse==1:
    num_iterations = np.round(interval_time / interval).astype(int) #ind the number of iterations for each interval as integer
    num_iterations[0]=num_iterations[0]+1


    current_timestamp=[] # counter for timestamp

    current_timestamp=0


    overall_start_time=time.time()
    
    if peltier==1:
        
        peltier_temp(pelt_temp_list[0], pelt_initial_heat_time)
        
        pelt_time=[] # counter for peltier time interval
        pelt_time=0
        ptc=0 # counter for pelt_Wait_schedule
        pelt_idx=0 # index for pelt_temp_list
        
        for interval_count in range(len(num_iterations)):
                for iteration in range(num_iterations[interval_count]):

                    if not (interval_count == 0 and iteration == 0):
                        current_timestamp = int(round(current_timestamp + interval[interval_count] * 60))
                        pelt_time = int(round(pelt_time + interval[interval_count]))
                        
                    index=0

                    start_time=time.time() # start of when imaging takes place

                    if pelt_idx==0:
                        peltier_temp(pelt_temp_list[pelt_idx], 1)
                        pelt_idx +=1
                    else:
                        if pelt_time >= pelt_wait_schedule[ptc]:
                            peltier_temp(pelt_temp_list[pelt_idx], 1)
                            pelt_idx +=1
                            pelt_time = 0
                            ptc += 1
                    

                    
                    for i in range(samples):

                        if enclosure_heating==1:

                            script_print_sample(f"sample_{sample_names[i]}_enclosure_temperature{temperature_enclosure}_peltier_temp:{pelt_temp_list[ptc]}_timestamp_{current_timestamp}_s", positions[index:], number[i])# gives the file a header and moves to positions and images sample
                            #print(current_timestamp)
                            index += number[i]
                        else: 
                            script_print_sample(f"sample_{sample_names[i]}_peltier_temp:{pelt_temp_list[ptc]}_timestamp_{current_timestamp}_s", positions[index:], number[i])# gives the file a header and moves to positions and images sample
                            #print(current_timestamp)
                            index += number[i]


                    end_time=time.time() # end of when imaging takes place

                    elapsed_time=(end_time-start_time) / 60 # time taken to image all samples in minutes

                    

                    if iteration<num_iterations[interval_count]-1:

                        interval_cycle=interval[interval_count]-elapsed_time # how long to wait till next imaging cycle
                        
                        # identify how long to sleep for
                        sleep_time = (timedelta(minutes=interval_cycle))

                        # Convert total seconds to hours, minutes, and seconds using map

                        rounded_seconds = round(sleep_time.total_seconds())
                        ht, mt, st = map(int, [
                            rounded_seconds // 3600,             # Hours
                            (rounded_seconds % 3600) // 60,      # Minutes
                            rounded_seconds % 60                # Seconds
                        ])
                        sleep_time_s = ht* 3600 + mt * 60 + st


                        # Countdown loop
                        for remaining in range(sleep_time_s, 0, -1):
                            mins, secs = divmod(remaining, 60)
                            hours, mins = divmod(mins, 60)
                            print(f"\rTime until next cycle: {hours:02}:{mins:02}:{secs:02}", end="")
                            time.sleep(1)
                            
                        print("\rTime until next cycle: 00:00:00")

                        print("\nImaging")
                    
                    else:
                        if interval[interval_count]!=interval[-1]:
                            interval_cycle=interval[interval_count+1]-elapsed_time # how long to wait till next imaging cycle
                            
                            # identify how long to sleep for
                            sleep_time = (timedelta(minutes=interval_cycle))

                            # Convert total seconds to hours, minutes, and seconds using map

                            rounded_seconds = round(sleep_time.total_seconds())
                            ht, mt, st = map(int, [
                                rounded_seconds // 3600,             # Hours
                                (rounded_seconds % 3600) // 60,      # Minutes
                                rounded_seconds % 60                # Seconds
                            ])
                            sleep_time_s = ht* 3600 + mt * 60 + st


                            # Countdown loop
                            for remaining in range(sleep_time_s, 0, -1):
                                mins, secs = divmod(remaining, 60)
                                hours, mins = divmod(mins, 60)
                                print(f"\rTime until next cycle: {hours:02}:{mins:02}:{secs:02}", end="")
                                time.sleep(1)
                                
                            print("\rTime until next cycle: 00:00:00")

                            print("\nImaging")



                        
                    

                        
    else:
        for interval_count in range(len(num_iterations)):
            for iteration in range(num_iterations[interval_count]):

                if interval_count == 0 and iteration==0:
                    current_timestamp=0 # counter for timestamp
                else:
                    current_timestamp=np.round((current_timestamp+(interval[interval_count]*60))).astype(int)
                
                index=0

                start_time=time.time() # start of when imaging takes place

                for i in range(samples):

                    if enclosure_heating==1:

                        script_print_sample(f"sample_{sample_names[i]}_enclosure_temperature{temperature_enclosure}_timestamp_{current_timestamp}_s_position_{index}", positions[index:], number[i])# gives the file a header and moves to positions and images sample
                        #print(current_timestamp)
                        index += number[i]
                    else:
                        script_print_sample(f"sample_{sample_names[i]}_timestamp_{current_timestamp}_s_position_{index}", positions[index:], number[i])# gives the file a header and moves to positions and images sample
                        #print(current_timestamp)
                        index += number[i]


                end_time=time.time() # end of when imaging takes place

                elapsed_time=(end_time-start_time) / 60 # time taken to image all samples in minutes

                if iteration<num_iterations[interval_count]-1:

                    interval_cycle=interval[interval_count]-elapsed_time # how long to wait till next imaging cycle
                    
                    # identify how long to sleep for
                    sleep_time = (timedelta(minutes=interval_cycle))

                    # Convert total seconds to hours, minutes, and seconds using map

                    rounded_seconds = round(sleep_time.total_seconds())
                    ht, mt, st = map(int, [
                        rounded_seconds // 3600,             # Hours
                        (rounded_seconds % 3600) // 60,      # Minutes
                        rounded_seconds % 60                # Seconds
                    ])
                    sleep_time_s = ht* 3600 + mt * 60 + st


                    # Countdown loop
                    for remaining in range(sleep_time_s, 0, -1):
                        mins, secs = divmod(remaining, 60)
                        hours, mins = divmod(mins, 60)
                        print(f"\rTime until next cycle: {hours:02}:{mins:02}:{secs:02}", end="")
                        time.sleep(1)
                        
                    print("\rTime until next cycle: 00:00:00")

                    print("\nImaging")
                else:
                        if interval[interval_count]!=interval[-1]:
                            interval_cycle=interval[interval_count+1]-elapsed_time # how long to wait till next imaging cycle
                            
                            # identify how long to sleep for
                            sleep_time = (timedelta(minutes=interval_cycle))

                            # Convert total seconds to hours, minutes, and seconds using map

                            rounded_seconds = round(sleep_time.total_seconds())
                            ht, mt, st = map(int, [
                                rounded_seconds // 3600,             # Hours
                                (rounded_seconds % 3600) // 60,      # Minutes
                                rounded_seconds % 60                # Seconds
                            ])
                            sleep_time_s = ht* 3600 + mt * 60 + st


                            # Countdown loop
                            for remaining in range(sleep_time_s, 0, -1):
                                mins, secs = divmod(remaining, 60)
                                hours, mins = divmod(mins, 60)
                                print(f"\rTime until next cycle: {hours:02}:{mins:02}:{secs:02}", end="")
                                time.sleep(1)
                                
                            print("\rTime until next cycle: 00:00:00")

                            print("\nImaging")





else: 
    if peltier ==1:
        peltier_temp(pelt_temp_list[0],1) #  only looks at 1st temperature in list as ramp is unnecessary
        index=0
        for i in range(samples):

            if enclosure_heating==1:

                script_print_sample(f"sample_{sample_names[i]}_enclosure_temperature{temperature_enclosure}_peltier_temp_{pelt_temp_list[0]}", positions[index:], number[i]) # gives the file a header and moves to positions and images sample
                index += number[i]
            else:
                script_print_sample(f"sample_{sample_names[i]}_peltier_temp_{pelt_temp_list[0]}", positions[index:], number[i]) # gives the file a header and moves to positions and images sample
                index += number[i]

    else:
        index=0
        for i in range(samples):

            if enclosure_heating==1:

                script_print_sample(f"sample_{sample_names[i]}_enclosure_temperature{temperature_enclosure}", positions[index:], number[i]) # gives the file a header and moves to positions and images sample
                index += number[i]
            else:
                script_print_sample(f"sample_{sample_names[i]}", positions[index:], number[i]) # gives the file a header and moves to positions and images sample
                index += number[i]
                

script_print_move_z_c_reset()


script_send_command(b"</temika>", False) # ends the temika script



############################################ Post Image processing #####################################################


### Create dictionary to find what each frame in each movie corresponds to
for sample_idx in range(len(sample_names)):
    # Build frame mapping
    frame_sequence = []
    frame_seq_id = 1

    if z_stack == 1:

        if z_c_order == 'zc':

            
                for num in range(1, number[sample_idx] + 1):
                    for zi in np.arange(-z_range,z_range +z_step,z_step):
                        for il in range(len(illumination)):
                            frame_sequence.append({
                                "frame": frame_seq_id,
                                "sample": sample_names[sample_idx],
                                "number": num,
                                "z": float(zi),
                                "illum_wavelength": illumination[il],
                                "illum_exposure_time": illum_expose[il],
                                "illum_pwr": laser_pwr[il]
                            })
                            frame_seq_id += 1

        elif z_c_order == 'cz':


                for num in range(1, number[sample_idx] + 1):
                    for il in range(len(illumination)):
                        for zi in np.arange(-z_range,z_range +z_step,z_step):                
                            frame_sequence.append({
                                "frame": frame_seq_id,
                                "sample": sample_names[sample_idx],
                                "number": num,
                                "z": float(zi),
                                "illum_wavelength": illumination[il],
                                "illum_exposure_time": illum_expose[il],
                                "illum_pwr": laser_pwr[il]
                            })

                            frame_seq_id += 1
    elif z_stack == 0:
        for num in range(1, number[sample_idx] + 1):
            for il in range(len(illumination)):              
                frame_sequence.append({
                    "frame": frame_seq_id,
                    "sample": sample_names[sample_idx],
                    "number": num,
                    "z": 1.0,
                    "illum_wavelength": illumination[il],
                    "illum_exposure_time": illum_expose[il],
                    "illum_pwr": laser_pwr[il]
                })

                frame_seq_id += 1

    csv_file = f"{main_folder}/frame_mapping_{sample_idx}.csv"
    with open(csv_file, mode="w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["frame", "sample", "number", "z", "illum_wavelength", "illum_exposure_time", "illum_pwr"])
        writer.writeheader()
        for f in frame_sequence:
            writer.writerow(f)

    print(f"Frame mapping export complete")

with open(f"{main_folder}/timestamp.csv", mode="w", newline="") as file:
    writer = csv.writer(file)
    for item in timestamp:
        writer.writerow([item])  # each item on a new row

print(f"timestamp export complete")


with open(f"{main_folder}/pos.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Timestamp", "x", "y", "z","a"])  # Add headers if needed
    writer.writerows(new_pos_save_vec)

print(f"z vec export complete")



######################## Variables to edit

movie_folder_path=movie_output_path


output_path=f'{main_folder}/tiffs'

os.makedirs(output_path, exist_ok=True)  # creates dircetory if missing

## Used this: illumination=[0x01,0x60,0x40,0x80] # change depending on which channel is being used

compression_type = "raw"  # or "tiff_lzw" for lossless compression or "raw" for no compresison



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

        if timelapse==1:
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
        else:
            for t in range(T):
                for c in range(C):
                    for z in range(Z):
                        planes.append(
                            Plane(
                                the_t=t,
                                the_c=c,
                                the_z=z
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

        if timelapse==1:

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
