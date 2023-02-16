# BAC Firing Detection Script
# Disclaimer: 'Interval' and 'Range' used as synonymes for the entire script

from pathlib import Path
from scipy.signal import argrelmax, argrelmin
import numpy as np
import matplotlib.pyplot as plt
import os


def calculate_highest_peak(chosen_array, parameter_highest_peak):
    """returns value of highest peak and its point in time (ms) - OPTIONAL FUNCTION"""
    if peak_orientation == "positive":
        peak_value = max
    else:   # peak_orientation == "negative":
        peak_value = min
    peak_index_x = np.where(chosen_array == peak_value(chosen_array[parameter_highest_peak]))[0][0]
    # Lines Below Optional And For Checking Purposes Only
    # print(f"Peak Amount: {len(chosen_array[parameter_highest_peak])}")
    # print("All Peaks:")
    # for i in range(0, len(chosen_array[parameter_highest_peak])):
    #     print(chosen_array[parameter_highest_peak[i]])
    print(f"Highest Peak: {peak_value(chosen_array[parameter_highest_peak])} At: {round(x[peak_index_x],3)}ms")
    return x[peak_index_x]

def calculate_ISI(chosen_array):
    """Na+/Soma ISI length and amount - REQUIRED FUNCTION"""
    ISI_list = []
    num_ISI = 0
    if len(chosen_array[parameter_peak]) == 1:
        print(f"Only One Peak At {round(x[parameter_peak[0]], 3)}ms\n" + "No ISI found")
    elif len(chosen_array[parameter_peak]) > 1:
        for num_peaks in range(0, len(chosen_array[parameter_peak]) - 1):
            ISI_value = x[parameter_peak[num_peaks + 1]] - x[parameter_peak[num_peaks]]
            num_ISI += 1
            print(f"ISI {num_ISI}) {round(ISI_value,3)}ms\n" + f"START: {round(x[parameter_peak[num_peaks]],3)}ms\n" + f"END: {round(x[parameter_peak[num_peaks + 1]],3)}ms")
    for num_x in range(0, len(chosen_array[parameter_peak])):
        ISI_list.append(x[parameter_peak[num_x]])
    return ISI_list

def calculate_Ca_peak_in_Na_ISI(x_ISI,x_peak_compare):
    """checks whether Ca2+/Apical Dendrite peak is located within Na+/Soma interval - OPTIONAL FUNCTION"""
    print("________________________________________________")
    if x_peak_compare <= x_ISI[0] or x_peak_compare >= x_ISI[len(x_ISI) - 1]:
        print("Ca2+ Channel/Apical Dendrite Peak is Outside ISI of Na+/Soma")
    else:
        for num in range(0,len(x_ISI) - 1):
            if x_ISI[num] <= x_peak_compare <= x_ISI[num+1]:
                print(("Ca2+ Channel/Apical Dendrite Peak is Between:\n") + (f"{round(x_ISI[num],3)}ms and {round(x_ISI[num+1],3)}ms"))

def interval_Ca(x_ISI):
    """returns Ca2+/Apical Dendrite interval and checks if it is in Na+/Soma interval - REQUIRED FUNCTION"""
    Ca_in_interval = False
    ms_filter = 200
    print("_______________________________________________________________________")
    print("Ca2+ Interval Scores:")
    print("---------------------")
# Threshold for Ca2+/Apical Dendrite Interval
    # Stopping Interval Determination When There Is Just A 'Straight Line'
    if calcium_min[0] == calcium_max[0]:
        return Ca_in_interval
    # Threshold Determination for NEGATIVE Orientated Graph
    elif peak_orientation == "negative":
        if min(calcium_min) < 0:
            threshold_Ca = BAC_array_1[4000] + (min(calcium_min) - BAC_array_1[4000]) * 0.25
        elif min(calcium_min) > 0:
            threshold_Ca = BAC_array_1[4000] - (BAC_array_1[4000] - min(calcium_min)) * 0.25
        elif min(calcium_max) > 0:
            threshold_Ca = BAC_array_1[4000] - (BAC_array_1[4000] - min(calcium_max)) * 0.25
        elif min(calcium_max) < 0:
            threshold_Ca = BAC_array_1[4000] + (min(calcium_max) - BAC_array_1[4000]) * 0.25
        print(f"Ca2+ Threshold to Determine Range: {threshold_Ca}")
        # Calculating Ca2+/Apical Dendrite Interval
        crossing_point = []
        for num in range(0, len(BAC_array_1) - 1):
            if num > ms_filter * 40:
                if BAC_array_1[num] > threshold_Ca >= BAC_array_1[num + 1]:
                    crossing_point.append(x[num + 1])
        for num in range(0, len(BAC_array_1) - 1):
            if num > ms_filter * 40:
                if BAC_array_1[num] <= threshold_Ca < BAC_array_1[num + 1]:
                    crossing_point.append(x[num])
        start_interval = crossing_point[0]
        end_interval = crossing_point[1]
        calcium_range = end_interval - start_interval
        print(f"Ca2+ Interval Start: {round(start_interval, 3)}ms")
        print(f"Ca2+ Interval End: {round(end_interval, 3)}ms")
        print(f"Distance Between 'Start' and 'End': {round(calcium_range, 3)}ms")
        print("_______________________________________________________________________")
        print("Comparing Na+/Soma Trace and Ca2+/Apical Dendrite Trace:")
        print("--------------------------------------------------------")
        print("Na+ Total Range (Including Tolerance Range):")
        pre_ISI = []
        post_ISI = []
        for num in range(0, len(BAC_array_2) - 1):
            if min_index_2[0] > num > ms_filter * 40:
                if BAC_array_2[num] < threshold_Ca <= BAC_array_2[num + 1]:
                    pre_ISI.append(x[num + 1])
        for num in range(0, len(BAC_array_2) - 1):
            if num > min_index_2[len(min_index_2) - 1]:
                if BAC_array_2[num] >= threshold_Ca > BAC_array_2[num + 1]:
                    post_ISI.append(x[num])
        print(f"START: {pre_ISI[0]}ms\nEND: {post_ISI[0]}ms")
        print(f"Distance Between 'Start' and 'End': {round(post_ISI[0] - pre_ISI[0], 3)}ms")
        # Checking When Ca2+ Channel/Apical Dendrite Interval Occurs
        if end_interval <= pre_ISI[0] or start_interval >= post_ISI[0]:
            print(
                "Ca2+ Channel/Apical Dendrite Trace is Completely Outside the Na+/Soma Range (ISI and Tolerance Range)")
            return Ca_in_interval
        else:
            # Checking START of Ca2+ Channel/Apical Dendrite Interval
            if pre_ISI[0] >= start_interval:
                print("Ca2+ Channel/Apical Dendrite Trace Starts Before Na+/Soma pre-ISI Tolerance Range at:\n" + (
                    f"{round(start_interval, 3)}ms"))
                Ca_in_interval = True
                # Starts Before pre-ISI
            if pre_ISI[0] <= start_interval <= x_ISI[0]:
                print("Ca2+ Channel/Apical Dendrite Trace Starts Inside Na+/Soma pre-ISI Tolerance Range at:\n" + (
                    f" {round(start_interval, 3)}ms"))
                Ca_in_interval = True
                # Starts Inside pre-ISI
            nISI = 1
            for num in range(0, len(x_ISI) - 1):
                if x_ISI[num] <= start_interval <= x_ISI[num + 1]:
                    nISI += num
                    print(("Ca2+/Apical Dendrite Trace Starts at:\n") + (
                        f"{round(start_interval, 3)}ms Within Na+/Soma ISI Number {nISI})"))
                    Ca_in_interval = True
                    # Starts Inside one of the ISIs
            if x_ISI[len(x_ISI) - 1] <= start_interval <= post_ISI[0]:
                print(("Ca2+/Apical Dendrite Trace Starts Inside Na+/Soma post-ISI Tolerance Range at:\n") + (
                    f"{round(start_interval, 3)}ms"))
                Ca_in_interval = True
                # Starts Inside post-ISI
            # Checking END of Ca2+ Channel/Apical Dendrite Interval
            nISI = 1
            for num in range(0, len(x_ISI) - 1):
                if x_ISI[num] <= end_interval <= x_ISI[num + 1]:
                    nISI += num
                    print(("Ca2+/Apical Dendrite Trace Ends at:\n") + (
                        f"{round(end_interval, 3)}ms Within Na+/Soma ISI Number {nISI})"))
                    Ca_in_interval = True
                    # Ends Inside one of the ISIs
            if x_ISI[len(x_ISI) - 1] <= end_interval <= post_ISI[0]:
                print(("Ca2+/Apical Dendrite Trace Ends Inside Na+/Soma post-ISI Tolerance Range at:\n") + (
                    f"{round(end_interval, 3)}ms"))
                Ca_in_interval = True
                # Ends Inside post-ISI
            if end_interval >= post_ISI[0]:
                print(("Ca2+/Apical Dendrite Trace Ends After Na+/Soma post-ISI Tolerance Range at:\n") + (
                    f"{round(end_interval, 3)}ms"))
                Ca_in_interval = True
                # Ends After post-ISI
            return Ca_in_interval
    # Threshold Determination for POSITIVE Orientated Graph
    elif peak_orientation == "positive":
        threshold_Ca = -65
        print(f"Ca2+ Threshold to Determine Interval: {threshold_Ca} mV")
        # Calculating Ca2+/Apical Dendrite Interval
        crossing_point = []
        for num in range(0, len(BAC_array_1) - 1):
            if num > ms_filter * 40:
                if BAC_array_1[num] < threshold_Ca <= BAC_array_1[num + 1]:
                    crossing_point.append(x[num + 1])
        for num in range(0, len(BAC_array_1) - 1):
            if num > ms_filter * 40:
                if BAC_array_1[num] >= threshold_Ca > BAC_array_1[num + 1]:
                    crossing_point.append(x[num])
        start_interval = crossing_point[0]
        end_interval = crossing_point[1]
        calcium_range = end_interval - start_interval
        print(f"Ca2+ Interval Start: {round(start_interval, 3)}ms")
        print(f"Ca2+ Interval End: {round(end_interval, 3)}ms")
        print(f"Distance Between 'Start' and 'End': {round(calcium_range, 3)}ms")
        print("_______________________________________________________________________")
        print("Comparing Na+/Soma Trace and Ca2+/Apical Dendrite Trace:")
        print("--------------------------------------------------------")
        print("Na+ Total Range (Including Tolerance Range):")
        pre_ISI = []
        post_ISI = []
        for num in range(0, len(BAC_array_2) - 1):
            if max_index_2[0] > num > ms_filter * 40:
                if BAC_array_2[num] < threshold_Ca <= BAC_array_2[num + 1]:
                    pre_ISI.append(x[num + 1])
        for num in range(0, len(BAC_array_2) - 1):
            if num > max_index_2[len(max_index_2) - 1]:
                if BAC_array_2[num] >= threshold_Ca > BAC_array_2[num + 1]:
                    post_ISI.append(x[num])
        print(f"START: {pre_ISI[0]}ms\nEND: {post_ISI[0]}ms")
        print(f"Distance Between 'Start' and 'End': {round(post_ISI[0]-pre_ISI[0],3)}ms")
        # Checking When Ca2+ Channel/Apical Dendrite Interval Occurs
        if end_interval <= pre_ISI[0] or start_interval >= post_ISI[0]:
            print("Ca2+ Channel/Apical Dendrite Trace is Completely Outside the Na+/Soma Range (ISI and Tolerance Range)")
            return Ca_in_interval
        else:
            # Checking START of Ca2+ Channel/Apical Dendrite Interval
            if pre_ISI[0] >= start_interval:
                print(("Ca2+ Channel/Apical Dendrite Trace Starts Before Na+/Soma pre-ISI Tolerance Range at:\n") + (f"{round(start_interval, 3)}ms"))
                Ca_in_interval = True
                # Starts Before pre-ISI
            if pre_ISI[0] <= start_interval <= x_ISI[0]:
                print(("Ca2+ Channel/Apical Dendrite Trace Starts Inside Na+/Soma pre-ISI Tolerance Range at:\n") + (f"{round(start_interval, 3)}ms"))
                Ca_in_interval = True
                # Starts Inside pre-ISI
            nISI = 1
            for num in range(0, len(x_ISI) - 1):
                if x_ISI[num] <= start_interval <= x_ISI[num + 1]:
                    nISI += num
                    print(("Ca2+/Apical Dendrite Trace Starts at:\n") + (f"{round(start_interval, 3)}ms Within Na+/Soma ISI Number {nISI})"))
                    Ca_in_interval = True
                    # Starts Inside one of the ISIs
            if x_ISI[len(x_ISI)-1] <= start_interval <= post_ISI[0]:
                print(("Ca2+ Channel/Apical Dendrite Trace Starts Inside Na+/Soma post-ISI Tolerance Range at:\n") + (f"{round(start_interval, 3)}ms"))
                Ca_in_interval = True
                # Starts Inside post-ISI
            # Checking END of Ca2+ Channel/Apical Dendrite Interval
            nISI = 1
            for num in range(0, len(x_ISI) - 1):
                if x_ISI[num] <= end_interval <= x_ISI[num + 1]:
                    nISI += num
                    print(("Ca2+ Channel/Apical Dendrite Trace Ends at:\n") + (f"{round(end_interval, 3)}ms Within ISI {nISI})"))
                    Ca_in_interval = True
                    # Ends Inside one of the ISIs
            if x_ISI[len(x_ISI) - 1] <= end_interval <= post_ISI[0]:
                print(("Ca2+ Channel/Apical Dendrite Trace Ends Inside Na+/Soma post-ISI Tolerance Range at:\n") + (f"{round(end_interval, 3)}ms"))
                Ca_in_interval = True
                # Ends Inside post-ISI
            if end_interval >= post_ISI[0]:
                print(("Ca2+/Apical Dendrite Trace Ends After Na+/Soma post-ISI Tolerance Range at:\n") + (f"{round(end_interval, 3)}ms"))
                Ca_in_interval = True
                # Ends After post-ISI
            return Ca_in_interval

def Na_Peaks():
    """Na+/Soma peak amount and values - REQUIRED FUNCTION"""
    if peak_orientation == "negative":
        threshold_Na = min(sodium_min) * 0.2
        Num_AP = []
        for i in range(0, len(sodium_min)):
            if sodium_min[i] < threshold_Na:
                Num_AP.append(sodium_min[i])
        AP_amount = len(Num_AP)
        return AP_amount
    if peak_orientation == "positive":
        threshold_Na = max(sodium_max) * 0.2
        Num_AP = []
        for i in range(0, len(sodium_max)):
            if sodium_max[i] > threshold_Na:
                Num_AP.append(sodium_max[i])
        AP_amount = len(Num_AP)
        return AP_amount

def points_on_graph(chosen_array, chosen_index_min, chosen_index_max):
    """shows peaks on the line - OPTIONAL FUNCTION"""
    for min_i in range(0, len(chosen_array[chosen_index_min])):
        plt.scatter(x[chosen_index_min[min_i]], chosen_array[chosen_index_min[min_i]],linewidth=0.3, s=25, c='r')
    for max_i in range(0, len(chosen_array[chosen_index_max])):
        plt.scatter(x[chosen_index_max[max_i]], chosen_array[chosen_index_max[max_i]],linewidth=0.3, s=25, c='g')

def threshold_Ca_peak():
    """sets threshold for Ca2+/Apical Dendrite - REQUIRED FUNCTION"""
    if peak_orientation == "positive":
        threshold_Ca = baseline2 + abs(baseline2 - max(sodium_max)) * 0.5
        if max(calcium_max) > threshold_Ca:
            return True
        else:
            return False
    elif peak_orientation == "negative":
        threshold_Ca = baseline2 - abs(baseline2 - min(sodium_min)) * 0.5
        if min(calcium_min) < threshold_Ca:
            return True
        else:
            return False


# Data Sets
BAC_list = ["Ca_HVA_apic2.txt", "NaTa_t_soma.txt", "apic2.txt", "soma.txt"]
label_list = ["Ca_HVA (Apical Dendrite)", "NaTa_t (Soma)", "Apical Dendrite", "Soma"]

# Asking Which Data Files Should Be Analyzed
time_condition = input(("Which Time Condition?\n") + ("Choose one of the following: 280 285 290 295 300 305 310 315 320\n"))
location = Path(__file__).parent
condition_folder_path = location / time_condition
folder_data_index = 0
subfolders = []
for f in os.scandir(condition_folder_path):
    if not f.name[0] == '.':
        subfolders.append(f.name)
        print(f"{folder_data_index}) {f.name}")
    folder_data_index += 1
setup = int(input(f"Choose Firing Setup! Type '0'-'{len(subfolders)-1}'\n"))
for files in range(len(BAC_list)):
    print(f"{files}) {BAC_list[files]}")
txt_file = int(input("Type '0' to Compare Current Traces, Or Type '2' to Compare Voltages Traces\n"))

# Finding Data File
pathBAC_1 = location / time_condition / subfolders[setup] / BAC_list[txt_file]
pathBAC_2 = location / time_condition / subfolders[setup] / BAC_list[txt_file + 1]

# Data From txt file Into List With Integers, And Later Into Array
# For FIRST Data Set
for num in open(pathBAC_1, 'r'):
    data_value = [eval(i) for i in num.split()]
BAC_array_1 = np.array(data_value)
# For SECOND Data Set
for num in open(pathBAC_2, 'r'):
    data_value = [eval(i) for i in num.split()]
BAC_array_2 = np.array(data_value)

# Defining x- And y-axis
x = []
x_interval = 0
for i in range(0, 20001):   # txt file with 20,000 data values into range of 0 to 500 for x-axis
    x_interval += 0.025
    x.append(x_interval)
y_1 = BAC_array_1
y_2 = BAC_array_2


# Minimum And Maximum Values Of FIRST Data Set
print("_______________________________________________________________________")
print(f"First Data Set: {BAC_list[txt_file]}")      # position 0 or 2 in BAC_list
print("-------------------------------------")
order_factor = 30                                   # to filter out peaks that are potentially interference factors
baseline1 = BAC_array_1[4000]
max_peak_y1 = BAC_array_1[np.argmax(BAC_array_1)]   # highest score of all MAXIMUM values
min_peak_y1 = BAC_array_1[np.argmin(BAC_array_1)]   # highest score of all MINIMUM values
print(f"Baseline: {baseline1}")
# Getting Indices Of MINIMUM Peak Value
if abs(baseline1 - max_peak_y1) > abs(baseline1 - min_peak_y1):
    peak_orientation = "positive"
    if min_peak_y1 > 0:
        threshold = min_peak_y1 * 0.4
    else:   # min_peak_y1 <= 0:
        threshold = min_peak_y1 * 1.6
    min_1 = BAC_array_1[argrelmin(BAC_array_1, order=order_factor)]
    min_index_1 = []
    for i in range(0, len(min_1)):
        if np.any(BAC_array_1[min_1[i] > threshold]):
            min_index_1.append(np.where(BAC_array_1 == min_1[i])[0][0])
    # print(f"Indices Of Min. Peaks: {min_index_1}")
else:
    peak_orientation = "negative"
    if min_peak_y1 > 0:
        threshold = min_peak_y1 * 1.6
    else:   # min_peak_y1 < 0:
        threshold = min_peak_y1 * 0.4
    min_1 = BAC_array_1[argrelmin(BAC_array_1, order=order_factor)]
    min_index_1 = []
    for i in range(0, len(min_1)):
        if np.any(BAC_array_1[min_1[i] < threshold]):
            min_index_1.append(np.where(BAC_array_1 == min_1[i])[0][0])
    # print(f"Indices Of Min. Peaks: {min_index_1}")
# Getting Values And Amount Of MINIMUM Peaks
calcium_min = []
if BAC_array_1[min_index_1].size == 0:
    calcium_min.append(BAC_array_1[np.argmin(BAC_array_1)])
    min_index_1.append(np.where(BAC_array_1 == BAC_array_1[np.argmin(BAC_array_1)])[0][0])
    print(f"Amount of Min. Peaks: {len(calcium_min)}")
    print("Minimum Peak Values:")
    print(calcium_min[0])
else:
    print(f"Amount of Min. Peaks: {len(BAC_array_1[min_index_1])}")
    print("Minimum Peak Values:")
    for i in range(0, len(BAC_array_1[min_index_1])):
        print(BAC_array_1[min_index_1[i]])
        calcium_min.append(BAC_array_1[min_index_1[i]])
# Getting Indices Of MAXIMUM Peak Value
if abs(baseline1 - max_peak_y1) > abs(baseline1 - min_peak_y1):
    peak_orientation = "positive"
    if max_peak_y1 > 0:
        threshold = max_peak_y1 * 0.4
    else:   # max_peak_y < 0
        threshold = max_peak_y1 * 1.6
    max_1 = BAC_array_1[argrelmax(BAC_array_1, order=order_factor)]
    max_index_1 = []
    for i in range(0, len(max_1)):
        if np.any(BAC_array_1[max_1[i] > threshold]):
            max_index_1.append(np.where(BAC_array_1 == max_1[i])[0][0])
    # print(f"Indices Of Max. Peaks: {max_index_1}")
else:
    peak_orientation = "negative"
    if max_peak_y1 > 0:
        threshold = max_peak_y1 * 1.6
    else:   # max_peak_y < 0
        threshold = max_peak_y1 * 0.4
    max_1 = BAC_array_1[argrelmax(BAC_array_1, order=order_factor)]
    max_index_1 = []
    for i in range(0, len(max_1)):
        if np.any(BAC_array_1[max_1[i] < threshold]):
            max_index_1.append(np.where(BAC_array_1 == max_1[i])[0][0])
    # print(f"Indices Of Max. Peaks: {max_index_1}")
# Getting Values And Amount Of MAXIMUM Peaks
calcium_max = []
if BAC_array_1[max_index_1].size == 0:
    calcium_max.append(BAC_array_1[np.argmax(BAC_array_1)])
    max_index_1.append(np.where(BAC_array_1 == BAC_array_1[np.argmax(BAC_array_1)])[0][0])
    print(f"Amount of Max. Peaks: {len(calcium_max)}")
    print("Maximum Peak Values:")
    print(calcium_max[0])
else:
    print(f"Amount of Max. Peaks: {len(BAC_array_1[max_index_1])}")
    print("Maximum Peak Values:")
    for i in range(0, len(BAC_array_1[max_index_1])):
        print(BAC_array_1[max_index_1[i]])
        calcium_max.append(BAC_array_1[max_index_1[i]])


# Minimum And Maximum Values Of SECOND Data Set
print("_______________________________________________________________________")
print(f"Second Data Set: {BAC_list[txt_file+1]}")   # position 1 or 3 in BAC_list
print("----------------------------------------")
order_factor = 30                                   # to filter out peaks that are potentially interference factors
baseline2 = BAC_array_2[4000]
max_peak_y2 = BAC_array_2[np.argmax(BAC_array_2)]   # highest score of all MAXIMUM values
min_peak_y2 = BAC_array_2[np.argmin(BAC_array_2)]   # highest score of all MINIMUM values
print(f"Baseline: {baseline2}")
# Getting Indices Of MINIMUM Peak Value
if abs(baseline2 - max_peak_y2) > abs(baseline2 - min_peak_y2):
    peak_orientation = "positive"
    if min_peak_y2 > 0:
        threshold = min_peak_y2 * 0.2
    else:
        threshold = min_peak_y2 * 1.8
    min_2 = BAC_array_2[argrelmin(BAC_array_2, order=order_factor)]
    min_index_2 = []
    for i in range(0, len(min_2)):
        if np.any(BAC_array_2[min_2[i] > threshold]):
            min_index_2.append(np.where(BAC_array_2 == min_2[i])[0][0])
    # print(f"Indices Of Min. Peaks: {min_index_2}")
else:
    peak_orientation = "negative"
    if min_peak_y2 > 0:
        threshold = min_peak_y2 * 1.8
    else:
        threshold = min_peak_y2 * 0.2
    min_2 = BAC_array_2[argrelmin(BAC_array_2, order=order_factor)]
    min_index_2 = []
    for i in range(0, len(min_2)):
        if np.any(BAC_array_2[min_2[i] < threshold]):
            min_index_2.append(np.where(BAC_array_2 == min_2[i])[0][0])
    # print(f"Indices Of Min. Peaks: {min_index_2}")
# Getting Values And Amount Of MINIMUM Peaks
sodium_min = []
if BAC_array_2[min_index_2].size == 0:
    calcium_min.append(BAC_array_2[np.argmin(BAC_array_2)])
    min_index_2.append(np.where(BAC_array_2 == BAC_array_2[np.argmin(BAC_array_2)])[0][0])
    print(f"Amount of Min. Peaks: {len(sodium_min)}")
    print("Minimum Peak Values:")
    print(sodium_min[0])
else:
    print(f"Amount of Min. Peaks: {len(BAC_array_2[min_index_2])}")
    print("Minimum Peak Values:")
    for i in range(0, len(BAC_array_2[min_index_2])):
        print(BAC_array_2[min_index_2[i]])
        sodium_min.append(BAC_array_2[min_index_2[i]])
# Getting Indices Of MAXIMUM Peak Value
if abs(baseline2 - max_peak_y2) > abs(baseline2 - min_peak_y2):
    peak_orientation = "positive"
    if max_peak_y2 > 0:
        threshold = max_peak_y2 * 0.2
    else:   # max_peak_y < 0
        threshold = max_peak_y2 * 1.8
    max_2 = BAC_array_2[argrelmax(BAC_array_2, order=order_factor)]
    max_index_2 = []
    for i in range(0, len(max_2)):
        if np.any(BAC_array_2[max_2[i] > threshold]):
            max_index_2.append(np.where(BAC_array_2 == max_2[i])[0][0])
    # print(f"Indices Of Max. Peaks: {max_index_2}")
else:
    peak_orientation = "negative"
    if max_peak_y2 > 0:
        threshold = max_peak_y2 * 1.8
    else:   # max_peak_y < 0
        threshold = max_peak_y2 * 0.2
    max_2 = BAC_array_2[argrelmax(BAC_array_2, order=order_factor)]
    max_index_2 = []
    for i in range(0, len(max_2)):
        if np.any(BAC_array_2[max_2[i] < threshold]):
            max_index_2.append(np.where(BAC_array_2 == max_2[i])[0][0])
    # print(f"Indices Of Max. Peaks: {max_index_2}")
# Getting Values And Amount Of MAXIMUM Peaks
sodium_max = []
if BAC_array_2[max_index_2].size == 0:
    calcium_max.append(BAC_array_2[np.argmax(BAC_array_2)])
    max_index_2.append(np.where(BAC_array_2 == BAC_array_2[np.argmax(BAC_array_2)])[0][0])
    print(f"Amount of Max. Peaks: {len(sodium_max)}")
    print("Maximum Peak Values:")
    print(sodium_max[0])
else:
    print(f"Amount of Max. Peaks: {len(BAC_array_2[max_index_2])}")
    print("Maximum Peak Values:")
    for i in range(0, len(BAC_array_2[max_index_2])):
        print(BAC_array_2[max_index_2[i]])
        sodium_max.append(BAC_array_2[max_index_2[i]])


# OPTIONAL: Highest Peak In Ca2+/Apical Dendrite Trace
print("_______________________________________________________________________")
print("Peak Score(s) for Ca2+ Channel/Apical Dendrite + Point in Time:")
print("---------------------------------------------------------------")
got_peaks = True
# Checks If There Are Any Peaks
if BAC_array_1[min_index_1].size == 0 and BAC_array_1[max_index_1].size == 0:
    if calcium_min[0] == calcium_max[0]:
        print(f"Peak Amount: 0")
        print("Graph Has No Amplitudes")
# Checking For Positive Orientated Peaks, If Not Then Calculating MINIMUM Peaks
elif BAC_array_1[max_index_1].size == 0:
    got_peaks = False
    print(f"No Maximum/Positive Peak")
    peak_orientation = "negative"
    highest_peak = calculate_highest_peak(chosen_array=BAC_array_1, parameter_highest_peak=min_index_1)
# Checking For Negative Orientated Peaks, If Not Then Calculating MAXIMUM Peaks
elif BAC_array_1[min_index_1].size == 0:
    got_peaks = False
    print(f"No Minimum/Negative Peak")
    peak_orientation = "positive"
    highest_peak = calculate_highest_peak(chosen_array=BAC_array_1, parameter_highest_peak=max_index_1)
# Determine Peak Orientation Before Executing Function 'calculate_highest_peak'
elif got_peaks:
    baseline = BAC_array_1[4000]
    max_peak_y = max(BAC_array_1[max_index_1])  # highest score of all MAXIMUM values
    min_peak_y = min(BAC_array_1[min_index_1])  # highest score of all MINIMUM values
    # Checks Which Peak Has Greater Distance From Baseline
    if abs(baseline-max_peak_y) > abs(baseline-min_peak_y):
        peak_orientation = "positive"
        index_highest_peak = max_index_1
        print("Peak Orientation: Positive")
    else:   # abs(baseline-max_peak_y) < abs(baseline-min_peak_y):
        peak_orientation = "negative"
        index_highest_peak = min_index_1
        print("Peak Orientation: Negative")
    highest_peak = calculate_highest_peak(chosen_array=BAC_array_1, parameter_highest_peak=index_highest_peak)


# Na+/Soma Inter-Spike-Interval (ISI) And Checking for BAC Firing
print("_______________________________________________________________________")
print("ISI Scores (Duration, Start, End):")
print("----------------------------------")
got_peaks = True
# Checks If There Are Any Peaks
if BAC_array_2[min_index_2].size == 0 and BAC_array_2[max_index_2].size == 0:
    print(f"Peak Amount: 0")
    print("Graph Has No Amplitudes")
# Checking For Positive Orientated Peaks, If Not Then Calculating MINIMUM Peaks
if BAC_array_2[max_index_2].size == 0:
    got_peaks = False
    print(f"No Maximum/Positive Peak")
    parameter_peak = min_index_2
    peak_orientation = "negative"
# Checking For Negative Orientated Peaks, If Not Then Calculating MAXIMUM Peaks
elif BAC_array_2[min_index_2].size == 0:
    got_peaks = False
    print(f"No Minimum/Negative Peak")
    parameter_peak = max_index_2
    peak_orientation = "positive"
# Determine Peak Orientation Before Checking For BAC Firing
elif got_peaks:
    baseline2 = BAC_array_2[4000]
    max_peak_y2 = max(BAC_array_2[max_index_2])     # highest score of all MAXIMUM values
    min_peak_y2 = min(BAC_array_2[min_index_2])     # highest score of all MINIMUM values
    if abs(baseline2-max_peak_y2) > abs(baseline2-min_peak_y2):
        peak_orientation = "positive"
        parameter_peak = max_index_2
        print("Peak Orientation: Positive")
        if len(BAC_array_2[parameter_peak]) == 1:
            print(f"Only One Peak At {round(x[parameter_peak[0]], 3)}ms\n" + "No ISI found")
    else:   # abs(baseline2-max_peak_y2) < abs(baseline2-min_peak_y2):
        peak_orientation = "negative"
        parameter_peak = min_index_2
        print("Peak Orientation: Negative")
        if len(BAC_array_2[parameter_peak]) == 1:
            print(f"Only One Peak At {round(x[parameter_peak[0]], 3)}ms\n" + "No ISI found")
    # Function Below Is Optional
    # calculate_Ca_peak_in_Na_ISI(x_ISI=calculate_ISI(BAC_array_2), x_peak_compare=highest_peak)

    # BAC Firing Criteria
    # Executing Function 'Na_Peaks()', 'interval_Ca' And 'threshold_Ca_peak()'
    if Na_Peaks() >= 2 and interval_Ca(x_ISI=calculate_ISI(BAC_array_2)) == True and threshold_Ca_peak():
        print("_______________________________________________________________________")
        print("Conclusion:")
        if peak_orientation == "negative":
            print("BAC Firing Is Likely")
            print("Check 'Voltage Trace' For Further Details.")
        else:
            print("BAC Firing")
    # Commands Below Are Optional And For Checking Purposes Only
    # elif Na_Peaks() >= 2:
    #     print("_______________________________________________________________________")
    #     print("Conclusion:")
    #     print("Probably Not BAC Firing")
    #     print("Only 'Na_Peaks() > 2' Criterion Met")
    # elif interval_Ca(x_ISI=calculate_ISI(BAC_array_2)) == True:
    #     print("_______________________________________________________________________")
    #     print("Conclusion:")
    #     print("Probably Not BAC Firing")
    #     print("Only 'interval_Ca(x_ISI=calculate_ISI(BAC_array_2)) == True' Criterion Met")
    else:
        print("_______________________________________________________________________")
        print("Conclusion:")
        if peak_orientation == "negative":
            print("Probably No BAC Firing")
            print("Check 'Voltage Trace' For Further Details.")
        else:
            print("No BAC Firing")
print("_______________________________________________________________________")


# PLOT
if txt_file == 0:
    plt.title(f"{time_condition}ms: Current Trace ({subfolders[setup]})")
else: # txt_file == 2
    plt.title(f"{time_condition}ms: Voltage Trace ({subfolders[setup]})")
plt.plot(x, y_1, color='blue', linestyle='solid', label=label_list[txt_file])
plt.plot(x, y_2, color='orange', linestyle='solid', label=label_list[txt_file+1])
plt.legend()
plt.xlabel("time (ms)")
if txt_file == 0:
    plt.ylabel("current (mA/cm²)")
else:   # txt_file == 2
    plt.ylabel("voltage (mV)")
plt.xticks(range(0,501,50))
plt.xlim([190, 510])
plt.plot(x,y_1)
points_on_graph(BAC_array_1, min_index_1, max_index_1)
points_on_graph(BAC_array_2, min_index_2, max_index_2)
plt.show()
