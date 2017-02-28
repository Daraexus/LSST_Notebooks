import matplotlib.pyplot as plt
import numpy as np


def get_frequency_dictionary(lcs, classification_dict, percentual=False):
    
    freq_dict = classification_dict.copy()
    for e in freq_dict.keys():
        freq_dict[e]=[]
    
    for lc in lcs:
        class_lc = lc.group_by("classification")
        indices = class_lc.groups.indices

        for key in freq_dict.keys():
            freq_dict[key].append(0)

        for i in range(1, len(indices)):

            i_start = indices[i-1] 
            i_end = indices[i]
            group = class_lc[i_start:i_end]
            if not percentual:
                freq_dict[group[-1]["classification"]][-1]=len(group)
            else:
                freq_dict[group[-1]["classification"]][-1]=float(len(group))/len(lc)*100
        
    return freq_dict



def plot_object_distribution(lcs, classification_dict):
    
    data = {}

    for lc in lcs:
        for point in lc:
            if data.has_key(point['classification']) == False:
                data[point['classification']] = 1
            else:
                data[point['classification']]+=1




    for key in classification_dict.keys():
        if data.has_key(key) == True:
            data[classification_dict[key]] = data.pop(key)

    
    
    plt.figure(figsize=(20,10))
    ind = range(0, len(data.values()))
    p = plt.bar(ind, data.values(),  align='center')
    plt.xticks(ind, data.keys())
    
    total = sum(data.values())

    for i, rect in enumerate(p):
        h = rect.get_height()
        perc = np.around(float(h)/total*100, 2)
        plt.text(rect.get_x()+ rect.get_width()/3, rect.get_height()+15, str(perc) +"% = "+  str(int(h)), fontsize=15)

    plt.show()
    
def plot_light_curve_with_tags(light_curve):
    
    plt.figure(figsize=(20,10))
    class_lc = light_curve.group_by("classification")
    indices = class_lc.groups.indices

    for i in range(1, len(indices)):
        i_start = indices[i-1] 
        i_end = indices[i]
        group = class_lc[i_start:i_end]
        marker = marker_dict[group[-1]["classification"]]
        label =  classification_dict[group[-1]["classification"]]
        plt.errorbar(group["mjd"], group["flux"],  yerr=group["flux_error"], fmt=marker, color='red', label=label, markersize=10)

    plt.legend(numpoints=1)
    #plt.ylim(ymin=0)

    plt.show()

def plot_proportion(x_data, y_data, x_data_label="x", y_data_label="y", percentual=False):
    
    
    
    plt.figure(figsize=(20,10))
    plt.scatter(x_data, y_data, 10)
    plt.xlabel(x_data_label)
    plt.ylabel(y_data_label)
    plt.xlim(xmin=-1)
    plt.ylim(ymin=-1)
    if percentual == False:
        max_n = max(max(x_data), max(y_data))+10
    else:
        max_n = 110
        plt.xticks(np.arange(0,101,10))
        plt.yticks(np.arange(0,101,10))
        plt.xlim(xmax=100)
        plt.ylim(ymax=100)
        
    plt.plot([0,max_n], [0,max_n], 'r--')
    
    for i in np.arange(0,max_n, 10):
        plt.plot([0,i], [i,0],  'r--')
    plt.grid(True)
    plt.show()
    
