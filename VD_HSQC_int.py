import numpy as np
import nmrglue as ng
import os, sys
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import tkinter as tk
from matplotlib.figure import Figure
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tqdm import tqdm
from tkinter import *
import importlib


# https://stackoverflow.com/questions/58367251/how-can-i-store-the-data-of-my-tkinter-entries-into-a-dataframe-to-later-export

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx, array[idx]

def Extract_NMR_Spectra_Bruker(
        path_nmr_data   =   'path_nmr_data',
        expno_data      =   'expno_data',
        procno_data     =   'procno_data',
        spec_lim        =   'spec_lim'
        ):
    path = os.path.join(path_nmr_data,str(expno_data),'pdata',str(procno_data))

    dic, data = ng.bruker.read_pdata(
        path,
        read_procs=True,
        read_acqus=False,
        scale_data = True,
        all_components=False
        )

    scaled_data = ng.fileio.bruker.scale_pdata(
        dic, 
        data, 
        reverse=False
        )
    if spec_lim == [0,0,0,0]:
        _data = scaled_data

    else:
        udic = ng.bruker.guess_udic(dic,data)
        uc_F1 = ng.fileiobase.uc_from_udic(udic, 0)
        uc_F2 = ng.fileiobase.uc_from_udic(udic, 1)
        ppm_scale_F2 = uc_F2.ppm_scale()
        ppm_scale_F1 = uc_F1.ppm_scale()
        idx_x0_F2, x0_F2 = find_nearest(ppm_scale_F2,spec_lim[0])
        idx_x1_F2, x1_F2 = find_nearest(ppm_scale_F2,spec_lim[1])
        idx_y0_F1, y0_F1 = find_nearest(ppm_scale_F1,spec_lim[2])
        idx_y1_F1, y1_F1 = find_nearest(ppm_scale_F1,spec_lim[3])
        _data = np.zeros((data.shape[0],data.shape[1]))
        data_ext = scaled_data[idx_y0_F1:idx_y1_F1,idx_x0_F2:idx_x1_F2]
        _data[idx_y0_F1:idx_y1_F1,idx_x0_F2:idx_x1_F2] = data_ext

    return _data, dic

def Peak_Picking_2D(
    data            =   'data', 
    threshold       =   'threshold',
    udic            =   'udic'
    ):
 
    peak_table = ng.peakpick.pick(
        data, 
        pthres=threshold, 
        algorithm='downward',
        )

    uc_F1 = ng.fileiobase.uc_from_udic(udic, 0)
    uc_F2 = ng.fileiobase.uc_from_udic(udic, 1)

    ppm_scale_F2 = uc_F2.ppm_scale()
    ppm_scale_F1 = uc_F1.ppm_scale()
    
    res = pd.DataFrame(columns=['ppm_F2_AXIS','ppm_F1_AXIS','pts_F2_AXIS','pts_F1_AXIS','peak_Amp'],index=np.arange(0,len(peak_table)-1))

    for l in range(len(peak_table)):
        res.loc[l,'pts_F2_AXIS'] = int(peak_table['X_AXIS'][l])+1
        res.loc[l,'pts_F1_AXIS'] = int(peak_table['Y_AXIS'][l])+1

        res.loc[l,'ppm_F2_AXIS'] = ppm_scale_F2[int(peak_table['X_AXIS'][l])]
        res.loc[l,'ppm_F1_AXIS'] = ppm_scale_F1[int(peak_table['Y_AXIS'][l])]

        peak_amplitudes = data[peak_table['Y_AXIS'][l].astype('int'),peak_table['X_AXIS'][l].astype('int')]
        res.loc[l,'peak_Amp'] = peak_amplitudes

    res = res.sort_values(by='pts_F1_AXIS', ascending=True)
    return res

def Plot_NMR_SubSpectra(data_All, pp_results, udic):
    uc_F1 = ng.fileiobase.uc_from_udic(udic, 0)
    uc_F2 = ng.fileiobase.uc_from_udic(udic, 1)

    ppm_scale_F2 = uc_F2.ppm_scale()
    ppm_scale_F1 = uc_F1.ppm_scale()
    
    F1_ppm_min = pp_results.ppm_F1_AXIS.min()-0.8
    F1_ppm_max = pp_results.ppm_F1_AXIS.max()+0.8

    F2_ppm_min = pp_results.ppm_F2_AXIS.min()-0.2
    F2_ppm_max = pp_results.ppm_F2_AXIS.max()+0.2

    idx_x0_F2, x0_F2 = find_nearest(ppm_scale_F2,F2_ppm_min)
    idx_x1_F2, x1_F2 = find_nearest(ppm_scale_F2,F2_ppm_max)

    idx_y0_F1, y0_F1 = find_nearest(ppm_scale_F1,F1_ppm_min)
    idx_y1_F1, y1_F1 = find_nearest(ppm_scale_F1,F1_ppm_max)
    
    SW_F1 = ppm_scale_F1.max()-ppm_scale_F1.min()
    TD_F1 = udic[0]['size']
    ppm_step_F1 = SW_F1/TD_F1

    fig = Figure(figsize=(10, 4), dpi=96)
    ax = fig.add_subplot(211)
    ax0 = fig.add_subplot(212)

    data_Pos_1k = data_All[2]
    data_Neg_1k = data_All[3]
    cl = 2e3*1.4**np.arange(30)

    cmap_list = ["viridis","plasma_r"]
    data_list = [data_Pos_1k,data_Neg_1k]

    data_Pos_16k = np.where(data_All[1]<0, 0, data_All[1])#data_All[4]
    data_Neg_16k = np.where(-data_All[1]<0, 0, -data_All[1])#data_All[5]
    data_list_1D = [data_Pos_16k,data_Neg_16k]
    _1D_location = int(pp_results.pts_F1_AXIS.mean())

    for k in range(len(data_list_1D)):
        data_1D = data_list_1D[k]
        _1D_Slice = data_1D[_1D_location,:]
        _1D_Slice = np.where(_1D_Slice<0, 0, _1D_Slice) 
        ax0.plot(
            ppm_scale_F2,
            _1D_Slice,
            color = matplotlib.cm.get_cmap(cmap_list[k])(0)
            )
    ax0.set_xlim(
        left = max(x1_F2,x0_F2),
        right = min(x1_F2,x0_F2)
    )


    for i in range(len(data_list)):
        data = data_list[i]
        subdata = data[idx_y1_F1:idx_y0_F1,idx_x1_F2:idx_x0_F2]
        ax.contour(
            subdata,
            cl,
            extent=(x1_F2,x0_F2,y1_F1,y0_F1),
            cmap = cmap_list[i]
            )
    ax.set_xlim(
        left = max(x1_F2,x0_F2),
        right = min(x1_F2,x0_F2)
    )
    ax.set_ylim(
        bottom = max(y1_F1,y0_F1),
        top = min(y1_F1,y0_F1)
    )
    Ylimits = ax.get_ylim()
    ax.set_xlabel(r'$^{1} H$ $(ppm)$')
    ax.set_ylabel(r'$^{13} C$ $(ppm)$')
    ax.xaxis.set_label_position('top')
    for i in range(len(pp_results)):
        ax.plot(
            pp_results.ppm_F2_AXIS.iloc[i],
            pp_results.ppm_F1_AXIS.iloc[i]-ppm_step_F1/2,
            marker='x',
            color='red'
            )
        ax.text(
            pp_results.ppm_F2_AXIS.iloc[i],
            Ylimits[0]-ppm_step_F1/2,
            'Peak '+str(i+1),
            rotation=90,
            verticalalignment='bottom'
            )

    return fig

def saveinfo(data, entries, entries_clustering):
    Peaks_Selection = []
    Clustering_Choice = []
    for i in range(len(entries)):
        Peaks_Selection.append(entries[i].get())
        Clustering_Choice.append(entries_clustering[i].get())

    data.append([Peaks_Selection,Clustering_Choice])
    
def opennewwindow(entries, entries_clustering, data, npeaks,pp_results, figure):

    # wdw.destroy()
    newwindow = tk.Tk()

    graph = FigureCanvasTkAgg(figure, master=newwindow)
    canvas = graph.get_tk_widget()
    canvas.grid(row=0, column=0,columnspan = 8,rowspan=1)

    data_cols_names = ['ppm_F2_AXIS','ppm_F1_AXIS','Spectra','Peaks_Selection','Cluster']

    for c in range(len(data_cols_names)):
        tk.Label(newwindow, text=str(data_cols_names[c]), ).grid(column=c+1, row=2)

    for i in range(npeaks):
        tk.Label(newwindow, text="Peak "+str(i+1), ).grid(column=0, row=i+3)

        en = tk.Entry(newwindow,justify = "center")
        en.insert(0, 'Yes')
        en.grid(column=len(data_cols_names)-1,row=i+3)
        entries.append(en)

        en_cluster = tk.Entry(newwindow,justify = "center")
        en_cluster.insert(0, 0)
        en_cluster.grid(column=len(data_cols_names),row=i+3,ipadx=5)
        entries_clustering.append(en_cluster)

        for c in range(len(data_cols_names)-2):
            col = data_cols_names[c]
            en_c = tk.Entry(newwindow,justify = "center")
            if col != 'Spectra':
                en_c.insert(0, round(pp_results.loc[i+1,col],2))
            if col == 'Spectra':
                en_c.insert(0, pp_results.loc[i+1,col])

            en_c.grid(column=c+1,row=i+3,sticky=tk.N+tk.S+tk.E+tk.W)
    

    tk.Button(newwindow, text="Save", command=lambda: saveinfo(data, entries, entries_clustering)).grid()
    tk.Button(newwindow, text="Exit", command=newwindow.destroy).grid()

    newwindow.mainloop()

def main_window(PeakPicking_data, figure):
    n_peak = len(PeakPicking_data)
    entries = []
    entries_clustering = []

    data = []

    
    opennewwindow(
        entries, 
        entries_clustering,
        data, 
        n_peak,
        PeakPicking_data, 
        figure
        )

    return data

def extract_1d(data, location, axis):
    """
    Extract a 1D slice from data along axis at location
    """
    s = [slice(v, v + 1) for v in location]
    s[axis] = slice(None, None)
    return np.atleast_1d(np.squeeze(data[tuple(s)]))

################################################################
# Test for missing librairies 
################################################################
if sys.argv[1] == "test":
    libfn = os.path.normpath(os.path.join(sys.argv[2],"vd_hsqc_lib.txt"))
    txt = ""
    modules = ["numpy","matplotlib","nmrglue","pandas","tkinter","tqdm"]
    for modname in modules:
        try:
            globals()[modname] = importlib.import_module(modname)
            txt += "\nLibrary installed : "+str(modname) 
        except ImportError as e:
            print(str(modname)+"is missing")
            txt += "\nLibrary missing : "+str(modname) 

    libf = open(libfn, 'w')
    libf.write(txt)
    libf.close()
    exit()
################################################################


topspin_dic = {
    'Path': str(sys.argv[1]),
    'Dummy_DataSet': int(sys.argv[5]), 
    'InPhase_DataSet': sys.argv[3], 
    'AntiPhase_DataSet': sys.argv[4], 
    'Data_Folder': sys.argv[2], 
    'Dummy_procno': int(sys.argv[6]),
    'Selected_window':[float(sys.argv[7]),float(sys.argv[8]),float(sys.argv[9]),float(sys.argv[10])]    
    }

################################################################
# Useful variables 
################################################################
Cleaning_threshold  =   0.33
Peaking_threshold   =   float(sys.argv[11])
y_Margin_Window     =   int(float(sys.argv[12]))
x_Margin_Window     =   int(float(sys.argv[13]))

################################################################
# Loading Data 
################################################################

# Created to later store all the data
Selected_Data_3D = []

# Path to export data back to bruker
path_out = os.path.join(topspin_dic['Path'],topspin_dic['Data_Folder'],str(topspin_dic['Dummy_DataSet']))

# Load, add, substract data and send data back to bruker
for i, (DataSet, ProcNo_In, ProcNo_Out) in enumerate(tqdm([
    (str(topspin_dic['InPhase_DataSet']),   1,      topspin_dic['Dummy_procno']     ),
    (str(topspin_dic['AntiPhase_DataSet']),   1,      topspin_dic['Dummy_procno']+1   ),
    (str(topspin_dic['AntiPhase_DataSet']),   1001,   topspin_dic['Dummy_procno']+2   ),
    (str(topspin_dic['AntiPhase_DataSet']),   1001,   topspin_dic['Dummy_procno']+3   ),
    ])):
    # print(ProcNo_Out)
    if i in [0,1,2,3]:
        data, dic_data = Extract_NMR_Spectra_Bruker(
                path_nmr_data   =   os.path.join(topspin_dic['Path'],topspin_dic['Data_Folder']),
                expno_data      =   DataSet,
                procno_data     =   ProcNo_In,
                spec_lim        =   topspin_dic['Selected_window']
                )
        dic_data_2_bruker = dic_data
        udic = ng.bruker.guess_udic(dic_data,data) # Need this universal directory for the peakc picking
        # plt.contour(
        #     data,
        #     2e3*1.4**np.arange(30)
        #     )
        # plt.show()
        # exit()
    if i in [0,1]: 
        # Read InPhase and AntiPhase dataset processed with 16k points
        Selected_Data_3D.append(data)

    if i == 2: 
        # Read AntiPhase dataset processed with 1k points and Set all negatives values to zero
        data = np.where(data<0, 0, data) 
        Selected_Data_3D.append(data)

    if i == 3:
        # Read AntiPhase dataset processed with 1k points and Set all positives values to zero
        data = np.where(-data<0, 0, -data) 
        Selected_Data_3D.append(data)
 
    # Export all datasets back to bruker 
    ng.fileio.bruker.write_pdata(
            dir = path_out,
            dic= dic_data,
            data= data,
            overwrite=True,
            pdata_folder=ProcNo_Out,       
            write_procs=False)
################################################################

################################################################
# Performing peak picking 
################################################################

# Running the peakc picking on both spectrum processed with 1k points in F2
appended_PeakPicking_Results = []
for k, (Description, data) in enumerate([
    ("Positive", Selected_Data_3D[2]),
    ("Negative", Selected_Data_3D[3]),
    ]):

    results = Peak_Picking_2D(
        data            =   data, 
        threshold       =   Peaking_threshold,
        udic            =   udic
        )
    results = results.assign(Spectra = str(Description)) 
    appended_PeakPicking_Results.append(results)

#DataFrame containing all peaks detected in the positive and the negative spectra.
appended_PeakPicking_Results = pd.concat(appended_PeakPicking_Results)
appended_PeakPicking_Results = appended_PeakPicking_Results.sort_values(by='pts_F1_AXIS', ascending=True)
appended_PeakPicking_Results.index = np.arange(1,len(appended_PeakPicking_Results)+1)

# Return a list of unique rows
unique_rows = appended_PeakPicking_Results['pts_F1_AXIS'].unique()

# 1st Cleaning procedure

Cleaned_PeakPicking_Results = pd.DataFrame()
for k in range(len(unique_rows)):
    PeakPicking_Results_Row = appended_PeakPicking_Results[appended_PeakPicking_Results['pts_F1_AXIS'] == unique_rows[k]]
    Sorted_PeakPicking_Results_Row = PeakPicking_Results_Row.sort_values(by='peak_Amp', ascending=False)

    n_peaks_F2 = len(PeakPicking_Results_Row)
    MeanAmp_PeakPicking_Results_Row = Sorted_PeakPicking_Results_Row.iloc[::1].loc[:,'peak_Amp'].mean()

    if n_peaks_F2 < 2:
        # remove rows for which only point in the F2 dimension is peaked
        pass

    if n_peaks_F2 == 2:
        Cleaned_PeakPicking_Results = Cleaned_PeakPicking_Results.append(PeakPicking_Results_Row)

    if n_peaks_F2 > 2:
        # Clean rows that contains more than 2 peaks assuming that all signals have to be at least more than 33% of average intensity of the most intense peaks of the row
        PeakPicking_Results_Row_selected = PeakPicking_Results_Row[PeakPicking_Results_Row.peak_Amp > MeanAmp_PeakPicking_Results_Row*Cleaning_threshold]
        Cleaned_PeakPicking_Results = Cleaned_PeakPicking_Results.append(PeakPicking_Results_Row_selected)

Cleaned_PeakPicking_Results = Cleaned_PeakPicking_Results.sort_values(by='pts_F1_AXIS', ascending=True)
################################################################

################################################################
# Matrix reconstruction
################################################################

# Data used for the shift
data_to_shift = Selected_Data_3D[0]

data_to_construct = np.zeros((Selected_Data_3D[0].shape[0],Selected_Data_3D[0].shape[1]))
cleaned_unique_rows = Cleaned_PeakPicking_Results['pts_F1_AXIS'].unique()



for i in range(len(cleaned_unique_rows)):
    pts_F1 = cleaned_unique_rows[i]
    PeakPicking_Results = Cleaned_PeakPicking_Results[Cleaned_PeakPicking_Results['pts_F1_AXIS'] == pts_F1]
    PeakPicking_Results = PeakPicking_Results.sort_values(by='ppm_F2_AXIS', ascending=True)
    PeakPicking_Results.index = np.arange(1,len(PeakPicking_Results)+1)
    if len(PeakPicking_Results) == 2 and len(np.unique(PeakPicking_Results['Spectra']))!= 1:
        pts_Offset = int(round((PeakPicking_Results['pts_F2_AXIS'].max()-PeakPicking_Results['pts_F2_AXIS'].min())/2))
        pts_Right_cmp = PeakPicking_Results.pts_F2_AXIS.max() #Right component
        data_to_shift_selected = data_to_shift[pts_F1-y_Margin_Window:pts_F1+y_Margin_Window,pts_Right_cmp-x_Margin_Window:pts_Right_cmp+x_Margin_Window]
        data_to_construct[pts_F1-y_Margin_Window:pts_F1+y_Margin_Window,pts_Right_cmp-x_Margin_Window-pts_Offset:pts_Right_cmp+x_Margin_Window-pts_Offset] += data_to_shift_selected

    else: 
        spec_fig = Plot_NMR_SubSpectra(Selected_Data_3D, PeakPicking_Results, udic)
        User_selection = main_window(PeakPicking_Results, spec_fig)
        if not User_selection: # If not choice is made by the user
            pass
        else:
            PeakPicking_Results['Peaks_Selection'] = User_selection[0][0]
            PeakPicking_Results['Cluster'] = User_selection[0][1]
            PeakPicking_Results = PeakPicking_Results[PeakPicking_Results.Peaks_Selection == "Yes"]

            n_cluster = PeakPicking_Results.Cluster.unique() # check the number of cluster

            for n in range(len(n_cluster)):

                Clustered_PeakPicking_Results = PeakPicking_Results[PeakPicking_Results.Cluster == n_cluster[n]]
                if len(Clustered_PeakPicking_Results) != 2 :
                    pass
                else:
                    pts_Offset = int(round((Clustered_PeakPicking_Results['pts_F2_AXIS'].max()-Clustered_PeakPicking_Results['pts_F2_AXIS'].min())/2))
                    pts_Right_cmp = Clustered_PeakPicking_Results.pts_F2_AXIS.max() #Right component
                    data_to_shift_selected = data_to_shift[pts_F1-y_Margin_Window:pts_F1+y_Margin_Window,pts_Right_cmp-x_Margin_Window:pts_Right_cmp+x_Margin_Window]
                    data_to_construct[pts_F1-y_Margin_Window:pts_F1+y_Margin_Window,pts_Right_cmp-x_Margin_Window-pts_Offset:pts_Right_cmp+x_Margin_Window-pts_Offset] += data_to_shift_selected
                    print(Clustered_PeakPicking_Results)
        


################################################################

################################################################
# Export Reconstructed matrix to bruker
################################################################

procno_shifted = topspin_dic['Dummy_procno']+4
ng.fileio.bruker.write_pdata(
        dir             =   path_out,
        dic             =   dic_data,
        data            =   data_to_construct,
        overwrite       =   True,
        pdata_folder    =   procno_shifted,       
        write_procs     =   False
        )
exit()

