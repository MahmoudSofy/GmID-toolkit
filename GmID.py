#Scripted by Mahmoud A. Sofy.
#Fayoum university, ECE class 2022.

###________________________________###

#V1.1 @ 2020
#General notes:
#lengths are mentioned in um, VDS and VBS are in mV, and rest of units are in SI.
#All functions and integrations were not set to handle wrong entries, so this toolkit is garbage in garbage out.
#Always check my github and linkedin since this script may be updated later.
#Your feedback will be welcomed!

###________________________________###

#Required libraries, make sure that you have all of them installed.
#You can use pip to check and install quickly numpy and matplotlib, rest are installed by default.

import numpy as np
import matplotlib.pyplot as plt
import math
import os
import re

###________________________________###

#First function: pre-conditioning the data, you may need to add more charachters depending on your simulator output.
#The same function can be done using macro on advanced editors like notepad++

def tabling(path):
    text = open(path)
    text = text.readlines()
    i = 0
    while i<len(text):
        if ((text[i].find('v')!=-1) or (text[i].find('l')!=-1)):
            del text[i]
            i = i-1
        i = i+1
    return np.genfromtxt(text)

###________________________________###

#Second function: creating data structure, only care about if your tables are not formatted as following:
#L(n):
    #VBS(m):
        #VGS(~|) | VDS(~-)
#You can use Create_DB() function to generate one file, you don't need to care about
#this function, except if your tables are not formatted as required.

def Create_DS(table, L_min = 0.2, L_step = 0.1, L_max = 3, VGS_min = 0, VGS_step = 0.025, VGS_max = 0.9, VDS_min = 0.1, VDS_step = 0.05, VDS_max = 0.9, VBS_min = 0, VBS_step = 0.1, VBS_max = 0.7):
    
    VGS_points = int((VGS_max - VGS_min)/VGS_step + 1)
    L_block = VGS_points * math.ceil((VBS_max - VBS_min)/VBS_step + 1)
    no_lengths = int(math.ceil(1 + (L_max-L_min)/L_step))
    no_VBS = int(math.ceil((VBS_max-VBS_min)/VBS_step))
    no_VDS = int(math.ceil((VDS_max-VDS_min)/VDS_step))
    
    #Initializing data blocks
    DS = {}
    flag = True
    for k in np.arange(0, no_lengths):
        for m in np.arange(0, no_VBS + 1):
            low = int(m*VGS_points + k*L_block)
            high = int(VGS_points + low)
            for j in np.arange(0, no_VDS + 1):
                atrI = 'VDS_' + str(int(1000*(VDS_min + j*VDS_step)))
                atrII = 'VBS_' + str(int(1000*(VBS_min + m*VBS_step)))
                if flag:
                    DS[atrI] = {}
                DS[atrI][atrII] = np.zeros([VGS_points,no_lengths])
            flag = False
    
    #creating data structure
    
    for k in np.arange(0, no_lengths):
        for m in np.arange(0, no_VBS + 1):
            low = int(m*VGS_points + k*L_block)
            high = int(VGS_points + low)
            for j in np.arange(0, no_VDS + 1):
                atrI = 'VDS_' + str(int(1000*(VDS_min + j*VDS_step)))
                atrII = 'VBS_' + str(int(1000*(VBS_min + m*VBS_step)))
                DS[atrI][atrII][:,int(k)] = table[low:high,int(j+1)]
    return DS
    
###________________________________###

#Third function: generating one file for each device containing all given parameters at the given corners
#Your files must be named as follow: devicetype + corner + parameter + .ext (e.g. pmosttgm.dat)
#Make sure you supply the function with the right sweeps.

def Create_DB(path, L_min = 0.2, L_step = 0.1, L_max = 3, VGS_min = 0, VGS_step = 0.025, VGS_max = 0.9, VDS_min = 0.1, VDS_step = 0.05, VDS_max = 0.9, VBS_min = 0, VBS_step = 0.1, VBS_max = 0.7):
    
    files = os.listdir(path)
    
    #Corners identification
    corner = []
    for i in np.arange(0, len(files)):
        if files[i].find('tt')!=-1:
            corner.append('TT')
        elif files[i].find('ff')!=-1:
            corner.append('FF')
        elif files[i].find('ss')!=-1:
            corner.append('SS')
    
    unique_corner = np.unique(np.array(corner))
    
    #Initializing data blocks
    DB = {}
    for i in np.arange(0, len(unique_corner)):
        DB[unique_corner[i]] = {}
        for k in np.arange(0, len(files)):
            DB[unique_corner[i]][files[k][6:-4]] = {}
        
    #DB creation
    for i in np.arange(0, len(files)):
        table = tabling(path+r'/'+files[i])
        DB[corner[i]][files[i][6:-4]] = Create_DS(table,  L_min = L_min, L_step = L_step, L_max = L_max, VGS_min = VGS_min, VGS_step = VGS_step, VGS_max = VGS_max, VDS_min = VDS_min, VDS_step = VDS_step, VDS_max = VDS_max, VBS_min = VBS_min, VBS_step = VBS_step, VBS_max = VBS_max)
    
    for i in np.arange(0, len(unique_corner)):
        DB[unique_corner[i]]['lengths'] = np.round(np.arange(L_min, L_max + L_step, L_step).reshape(1, -1), 2)
        DB[unique_corner[i]]['vgs'] = np.arange(VGS_min, VGS_max + VGS_step, VGS_step)
    return DB
    
###________________________________###

#Fourth function: generate the required charts
#Device width used here is 2u/finger * 5 fingers, change these values if you have used different values.

def charting(data_file, lengths, x_range, corner = 'TT', VDS = 600, VBS = 0, x = 'gm/id', y = ['gm/gds', 'id/w', 'id/gds', 'vgs', 'cdd/cgg', 'vdsat', 'gm/(cdd+cgg)']):
    
    width = 2e-6
    no_fingers = 5
    w = width * no_fingers
    
    #identifying the required parameters to be loaded
    parameters = []
    for i in np.arange(0, len(y)):
        parameters = list(set(parameters) | set([i for i in re.split(r'[/,(,),+,-,*]', y[i]) if i]))
    
    xparameters = set([i for i in re.split(r'[/,(,),+,-,*]', x) if i])
    for new in xparameters:
        if new not in parameters:
            parameters.append(new)
    try:
        parameters.remove('w')
    except:
        pass
    
    #identifying the required lengths
    temp1, colomns, temp2 = np.intersect1d(data_file[corner]['lengths'], lengths, return_indices=True)
    
    #Fetching the required tables
    vds = 'VDS_' + str(VDS)
    vbs = 'VBS_' + str(VBS)
    for i in np.arange(0, len(parameters)):
        try:
            exec(
                 parameters[i] + '=' +'data_file[corner][parameters[i]][vds][vbs][:, colomns]'
                )
        except:
            exec(
                 parameters[i] + '=' +'data_file[corner][parameters[i]]'
                 )
            exec(
                parameters[i] + '= np.outer(np.ones(len(colomns)), '+parameters[i]+').transpose()'
                )
            
    #pre-ploting conditioning
    subplots = len(y)
    subplots_rows = math.floor(math.sqrt(subplots))
    subplots_colomns = math.ceil(math.sqrt(subplots))
    if (subplots_colomns*subplots_rows) < subplots:
        subplots_colomns = subplots_colomns +1
    legends = []
    for i in lengths:
        legends.append(str(i)+'um')
        
    #Plotting
    ax = [None] * subplots
    xaxis = abs(eval(x))
    pts = np.where((xaxis <= x_range[1]) & (xaxis >= x_range[0]))
    fig = plt.figure()
    for i in np.arange(0, subplots):
        ax[i] = fig.add_subplot(subplots_rows, subplots_colomns, i+1)
        yaxis = abs(eval(y[i]))
        ax[i].plot(xaxis,yaxis, linewidth=0.8)
        ax[i].set_xlim(x_range)
        ymin = np.min(yaxis[pts])
        ymax = np.max(yaxis[pts])
        ax[i].set_ylim([ymin, ymax])
        if (ymax-ymin)>1e4:
            plt.yscale('log')
        ax[i].set(xlabel='$'+x+'$', ylabel='$'+y[i]+'$')
        ax[i].yaxis.set_label_coords(1, 0.5)
        ax[i].xaxis.set_label_coords(0.5, 1)
        ax[i].tick_params(axis='both', which='major', labelsize=8)
        plt.minorticks_on()
        ax[i].grid(which='major', linestyle='-', linewidth=0.5)
        ax[i].grid(which='minor', linestyle=':', linewidth=0.4)
        
    #Post plotting conditioning
    fig.legend(legends, loc='upper right')
    fig.suptitle('Vds = '+str(VDS)+'mV, Vbs = '+str(VBS)+'mV, @'+corner, fontsize=12)
    plt.figtext(0.005,0.005,'Scripted by Mahmoud A. Sofy', fontsize='xx-small')
    plt.show(block=True)

###________________________________###

#Fifth function: visualize the effect of Vds on the different parameters
#If the previous functions worked properly this will work properly also

def charting3D(data_file, length, x='gm/id', y = ['gm/gds', 'id/w', 'id/gds', 'vgs', 'cdd/cgg', 'vdsat', 'gm/(cdd+cgg)'], corner = 'TT', VBS = 0):
    
    width = 2e-6
    no_fingers = 5
    w = width * no_fingers
    
    #identifying the required parameters to be loaded
    parameters = []
    for i in np.arange(0, len(y)):
        parameters = list(set(parameters) | set([i for i in re.split(r'[/,(,),+,-,*]', y[i]) if i]))
        
    xparameters = re.split(r'[/,(,),+,-,*]', x)
    for new in xparameters:
        if new not in parameters:
            parameters.append(new)
    try:
        parameters.remove('w')
    except:
        pass

    #identifying the required length
    temp1, colomn, temp2 = np.intersect1d(data_file[corner]['lengths'], length, return_indices=True)
    
    #fetching the required tables
    vdss = list(data_file['TT']['id'])
    VBS = 'VBS_' + str(VBS)
    vgs_points = len(data_file[corner]['vgs'])
    for j in np.arange(0,len(parameters)):
        exec(
            parameters[j] + '= np.zeros([vgs_points,1])'
            )
        for i in np.arange(0,len(vdss)):
            try:
                exec(
                    parameters[j] + '= np.column_stack((eval(parameters[j]), data_file[corner][parameters[j]][vdss[i]][VBS][:,colomn]))'
                    )
            except:
                exec(
                    parameters[j] + '= np.column_stack((eval(parameters[j]), data_file[corner][parameters[j]]))'
                    )
        exec(
            parameters[j] + '= np.delete(eval(parameters[j]), 0, 1)'
            )
    
    #pre-ploting conditioning
    subplots = len(y)
    subplots_rows = math.floor(math.sqrt(subplots))
    subplots_colomns = math.ceil(math.sqrt(subplots))
    if (subplots_colomns*subplots_rows) < subplots:
        subplots_colomns = subplots_colomns + 1
    legends = []
    for i in vdss:
        legends.append(i.replace('_', '=')+'mV')
        
    #plotting
    ax = [None] * len(y)
    xaxis = abs(eval(x))
    fig = plt.figure()
    for i in np.arange(0, subplots):
        ax[i] = fig.add_subplot(subplots_rows, subplots_colomns, i+1, projection='3d')
        yaxis = abs(eval(y[i]))
        for j in np.arange(0, len(vdss)):
            z = int(vdss[j].split('_')[1])
            ax[i].plot(xaxis[:,j], yaxis[:,j], z, zdir='y')
        ax[i].set(xlabel='$'+x+'$', zlabel='$'+y[i]+'$', ylabel='VDS')
 
    #post plotting
    fig.legend(legends, loc='upper right', fontsize=7)
    fig.suptitle('L = '+str(length) +'um, Vbs = '+str(VBS.split('_')[1])+'mV, @'+corner, fontsize=12)
    plt.figtext(0.005,0.005,'Scripted by Mahmoud A. Sofy', fontsize='xx-small')
    plt.show(block=True)

###________________________________###

#INTEGRATION

logo = '''

 ██████╗ ███╗   ███╗██╗██████╗        ██╗   ██╗ ██╗   ██╗
██╔════╝ ████╗ ████║██║██╔══██╗       ██║   ██║███║  ███║
██║  ███╗██╔████╔██║██║██║  ██║       ██║   ██║╚██║  ╚██║
██║   ██║██║╚██╔╝██║██║██║  ██║       ╚██╗ ██╔╝ ██║   ██║
╚██████╔╝██║ ╚═╝ ██║██║██████╔╝███████╗╚████╔╝  ██║██╗██║
 ╚═════╝ ╚═╝     ╚═╝╚═╝╚═════╝ ╚══════╝ ╚═══╝   ╚═╝╚═╝╚═╝
                                                         
████████╗ ██████╗  ██████╗ ██╗     ██╗  ██╗██╗████████╗  
╚══██╔══╝██╔═══██╗██╔═══██╗██║     ██║ ██╔╝██║╚══██╔══╝  
   ██║   ██║   ██║██║   ██║██║     █████╔╝ ██║   ██║     
   ██║   ██║   ██║██║   ██║██║     ██╔═██╗ ██║   ██║     
   ██║   ╚██████╔╝╚██████╔╝███████╗██║  ██╗██║   ██║     
   ╚═╝    ╚═════╝  ╚═════╝ ╚══════╝╚═╝  ╚═╝╚═╝   ╚═╝     
                                                                                                     
'''
tools = '''
____________| Data files |____________
[1]Create data file

_____________| Charting |_____________
[2]Generate fundamental GmID charts
                                    
_________| Advanced Charting |________
[3]Advanced charts      [4]3D charts

______________| EXTRAS |______________
[5]Help      [6]Contact me
'''
sep = '''
______________________________________
'''
cp = os.getcwd()
info = '''
Mahmoud A. Sofy
Fayoum university, ECE class 2022
______________________________________
Email: ma3644@fayoum.edu.eg
Linkedin: /in/mahmoud-sofy
'''
Help = '''
Welcome to V1.0:
    
Overview:
   This script uses data tables to plot gm/id charts mainly
   but you can use advanced charting to plot any parameters
   and figure out how your trade off is.
   
[1]Create data file:
    First we need to have data, you can use any simulator and
    any model, but you have to care about the output data format
    the format which this script supports, for now, is:
        length(~|)
            vbs(~|)
                vgs(~|) vds(~-)
    for cadence virtuoso you can use adexl to export ocean script, 
    take caring of saving the required parameters, and the sweeps
    to be done in the same sequence, so you set the first variable
    to sweep is length and so on.
    Then run the ocean script using the log window and then export 
    all parameters you need, each in individual file.
    
    Then we need to convert this file to a table which this script
    can use, so we use this tool "Create data file" to generate this
    table and choose a suitable name to be easy to know the device
    type and technology later.
    
[2]Generate fundamental GmID charts:
    This tool will generate the basic charts you need, with minimum
    effort, all you need is to determine the basics.
    
[3]Advanced charts:
    This tool will let you to choose more than the basics, you can
    choose x-axis and y-axes so that you can plot more charts than
    the basics, and also you can choose another x-axis so you can
    visualize your trade offs better.
    
[4]3D charts:
    This tool will unleash your insight, you can visualize your
    device at different VDS, and how these curves vary, you can
    plot anything versus anything, so you can imagine how your
    device will behave in different cases.
    
**contact me:
    Your feedback will be welcomed, and you can ask me if you
    need.
'''
print(logo)

while True:
    print(tools)
    tool = input('>>> ')
    if tool == '1':
        file = Create_DB(input('Path of the parameters files>>> '), L_min = float(input('Minimum length>>> ') or 0.2), L_step = float(input('Length step>>> ') or 0.1), L_max = float(input('Maximum length>>> ') or 3), VGS_min = float(input('Minimum VGS>>> ') or 0), VGS_step = float(input('VGS step>>> ') or 0.025), VGS_max = float(input('Maximum VGS>>> ') or 0.9), VDS_min = float(input('Minimum VDS>>> ') or 0.1), VDS_step = float(input('VDS step>>> ') or 0.05), VDS_max = float(input('Maximum VDS>>> ') or 0.9), VBS_min = float(input('Minimum VBS>>> ') or 0), VBS_step = float(input('VBS step>>> ') or 0.1), VBS_max = float(input('Maximum VBS>>> ') or 0.7))
        path = input('Save the file as>>> ') + '.npy'
        np.save(path, file)

    elif tool == '2':
        available_files = os.listdir(cp)
        print('Available files: \n', available_files)
        filein = input('Which data file do you need>>> ')
        data_file = np.load(filein, allow_pickle=True)[()]
        lengthsin = eval(input('Enter the required [lengths]>>> ')) or [0.2, 0.4, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3]
        limitsin = eval(input('Enter the [limits]>>> ')) or [5, 20]
        cornerin = input('Enter the corner>>> ') or 'TT'
        VDSin = int(input('Enter the VDS>>> ') or 600)
        VBSin = int(input('Enter the VBS>>> ') or 0)
        xin = 'gm/id'
        yin = ['gm/gds', 'id/w', 'id/gds', 'vgs', 'cdd/cgg', 'vdsat', 'gm/(cdd+cgg)']
        charting(data_file, lengthsin, limitsin, corner = cornerin, VDS = VDSin, VBS = VBSin, x = xin, y = yin)
        
    elif tool == '3':
        available_files = os.listdir(cp)
        print('Available files: \n', available_files)
        filein = input('Which data file do you need>>> ')
        data_file = np.load(filein, allow_pickle=True)[()]
        lengthsin = eval(input('Enter the required [lengths]>>> '))
        limitsin = eval(input('Enter the [limits]>>> ')) or [5, 20]
        cornerin = input('Enter the corner>>> ') or 'TT'
        VDSin = int(input('Enter the VDS>>> ') or 600)
        VBSin = int(input('Enter the VBS>>> ') or 0)
        xin = input('Enter the x-axis>>> ') or 'gm/id'
        yin = eval(input('Enter the [y-axes]>>> ') or '''['gm/gds', 'id/w', 'id/gds', 'vgs', 'cdd/cgg', 'vdsat', 'gm/(cdd+cgg)']''')
        charting(data_file, lengthsin, limitsin, corner = cornerin, VDS = VDSin, VBS = VBSin, x = xin, y = yin)
       
    elif tool == '4':
        available_files = os.listdir(cp)
        print('Available files: \n', available_files)
        filein = input('Which data file do you need>>> ')
        data_file = np.load(filein, allow_pickle=True)[()]
        lengthsin = eval(input('Enter the required length>>> '))
        cornerin = input('Enter the corner>>> ') or 'TT'
        VBSin = int(input('Enter the VBS>>> ') or 0)
        xin = input('Enter the x-axis>>> ') or 'gm/id'
        yin = eval(input('Enter the [y-axes]>>> ') or '''['gm/gds', 'id/w', 'id/gds', 'vgs', 'cdd/cgg', 'vdsat', 'gm/(cdd+cgg)']''')
        charting3D(data_file, lengthsin, corner = cornerin, VBS = VBSin, x = xin, y = yin)
    
    elif tool == '5':
        print(Help)
        
    elif tool == '6':
        print(info)
    else:
        print('Enter a valid choice')
    print(sep)