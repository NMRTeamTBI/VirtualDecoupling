#########
# Virtual Decoupling
#########

#########
# This script reads a pair of In-phase/Anti-phase spectra
#########

###############
# Changes log #
###############
#		v0.1 (26/01/20)
#			-   initial script
#       v0.2 (01/02/20)
#           -   Perform Peak picking on the inphase dataset 
#       v0.3 (05/02/20)
#           -   Perform Peak picking on the antiphase dataset 
#               (positive and negative)
#           -   Implementation of the user interface to select peaks to consider
#               and cluster them (Visually select which right and left component)
#       v0.4 (22/03/20)
#           -   Implementation of the choice for python path
###############

##########
# Script #
##########

import os, sys, subprocess
from shutil import copyfile, copytree, rmtree
from datetime import *
import importlib

def initialize(options):
	try:
		options = readOpt(options)
	except:
		writeOpt(options)
	return(options)

def readOpt(options):
	optfn = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), "vd_hsqc_path.txt"))
	optf = open(optfn, 'r')
	for l in optf.readlines():
		li = l.strip("\n").split("\t")
		if li[0] in options.keys():
		  options[li[0]] = li[1]
	optf.close()
	options = parseOpt(options)
	return(options)

def writeOpt(options):
  optfn = os.path.normpath(os.path.join(os.path.dirname(sys.argv[0]), "vd_hsqc_path.txt"))
  optt = "\n".join([k + "\t" + str(v) for k, v in options.items()])
  optf = open(optfn, 'w')
  optf.write(optt)
  optf.close()

def parseOpt(options):
	if options["python_path"] != "python" and not os.path.isfile(options["python_path"]):
	  ERRMSG(message = "Python path not found.", title="Error", details=None, modal=1)
	  EXIT()
	return(options)

def modifyOpt(options):
    optold = initialize(options)
    diagOpt = INPUT_DIALOG("Virtual Decoupling", "Modify path options.", ["python path", "script path"], [optold["python_path"],optold["script_path"]], ["",""], ["1","1"])
    if diagOpt != None:
        kopt = ["python_path","script_path"]
        for o,v in enumerate(kopt):
            optold[v] = diagOpt[o]
        opt_parsed = parseOpt(optold)
        writeOpt(opt_parsed)
    EXIT()

def copy_all_pdata_files(src,dest):
    src_files = os.listdir(src)
    for file_name in src_files:
        src_file_name = os.path.normpath(os.path.join(src, file_name))
        dest_file_name = os.path.normpath(os.path.join(dest, file_name))
        if os.path.isfile(src_file_name):
            copyfile(src_file_name, dest_file_name)

today = date.today()
day = today.strftime("%d/%m/%Y")
now = datetime.now()
current_time = now.strftime("%H:%M:%S")

#######################################################
# default parameters
opt_ini = {"python_path":"python", "script_path":str(os.path.dirname(sys.argv[0]))}
# load user-defined parameters (default values are used if none declared)
options = initialize(opt_ini) 
python_env = options["python_path"]

#######################################################
# check python path and the path for the analysis script
if "--opt" in sys.argv:
    modifyOpt(options)
    EXIT()

Analysis_Script_Path = os.path.normpath(os.path.join(options["script_path"],"VD_HSQC_int.py"))

#######################################################
# check if all python dependencies are available for the analysis script
if "--test" in sys.argv:
#     cmd = "python /opt/topspin4.0.8/exp/stan/nmr/py/user/VD_HSQC_int.py %s %s"  % ('test', options["script_path"])
    cmd = str(os.path.normpath(os.path(options["python_path"])))+" "+str(Analysis_Script_Path)+" %s %s"  % ('test', options["script_path"])
    subprocess.Popen(cmd, shell=True).wait()
    EXIT()

#######################################################
# check if multiple display is active
if SELECTED_WINDOW().isMultipleDisplayActive():
    ERRMSG(message = "Please exit multiple display before running the script.", title="Error", details=None, modal=1)
    EXIT()

# get current dataset
current_dataset = CURDATA()

#######################################################
# ask for imputs
inputs = INPUT_DIALOG("Virtual Decoupling",\
"Please make sure that you New ExpNo does not correspond to an experimental dataset \n",
# "\n",
# "New ExpNo : ExpNo that will contain the resulting spectra ",
["Data Directory","InPhase ExpNo ", "AntiPhase ExpNo", "New ExpNo ", "First ProcNo ", "Chemical shift(s) in ppm","Comment (Title)", "nc_proc","Detection threshold","Points Shifting"], [str(current_dataset[0]),"3", "5", "776", "1","3.45;3.1;39.5;37.5","","-4", "0.1e4","10;40"],\
["","","","ExpNo that will contain the resulting spectra","First ProcNo. The following ones will be incremented automatically","max(F2);min(F2);max(F1);min(F1)  // Leave zeros for full spectrum","Any comment(s) that should be added to the title","nc proc parameters for visualization (default)","Threshold for peak picking (default)","Number of points for box selection (default)"], ["1", "1", "1", "1","1"])

if len(inputs[0]) == 0 or len(inputs[1]) == 0:
    ERRMSG(message = "Please provide the expnos for the InPhase and AntiPhase experiments ", title="Error", details=None, modal=1)
else:
    expno_InPhase = int(inputs[1])
    expno_AntiPhase = int(inputs[2])
    dummy_expno = int(inputs[3])
    dummy_procno = int(inputs[4])
    Comment_Title = str(inputs[6])
    nc_proc = str(inputs[7])
    pp_th = str(inputs[8])
    Box_pts = [float(i) for i in inputs[9].split(";")]
    Wdw_ppm = [float(i) for i in inputs[5].split(";")]

#######################################################
# Check that both InPhase and AntiPhase spectra have the same dimennsions (STSR & STSI should be the same in both)
RE([str(current_dataset[0]), str(expno_AntiPhase), "1", str(current_dataset[3])+"/"])
STSR_AntiPhase_F2 = GETPAR("STSR",axis=0)
STSI_AntiPhase_F2 = GETPAR("STSI",axis=0)
STSR_AntiPhase_F1 = GETPAR("STSR",axis=1)
STSI_AntiPhase_F1 = GETPAR("STSI",axis=1)

RE([str(current_dataset[0]), str(expno_InPhase), "1", str(current_dataset[3])+"/"])
STSR_InPhase_F2 = GETPAR("STSR",axis=0)
STSI_InPhase_F2 = GETPAR("STSI",axis=0)
STSR_InPhase_F1 = GETPAR("STSR",axis=1)
STSI_InPhase_F1 = GETPAR("STSI",axis=1)

if (STSR_AntiPhase_F2 != STSR_InPhase_F2) or (STSI_AntiPhase_F2 != STSI_InPhase_F2) or (STSR_AntiPhase_F1 != STSR_InPhase_F1) or (STSI_AntiPhase_F1 != STSI_InPhase_F1):
    ERRMSG(message = "Both the InPhase and AntiPhase spectra should have the same window \n (STSR & STSI should be equals both datasets)", title="Error", details=None, modal=1)

#######################################################
# Set nc_proc identically in both InPhase
RE([str(current_dataset[0]), str(expno_InPhase), "1", str(current_dataset[3])+"/"])
XCMD('xfb nc_proc '+str(nc_proc))

RE([str(current_dataset[0]), str(expno_AntiPhase), "1", str(current_dataset[3])+"/"])
XCMD('xfb nc_proc '+str(nc_proc))

#######################################################
# Create a dummy directory to store the additional spectra 
path_dummy_expno = os.path.normpath(os.path.join(current_dataset[3], str(inputs[0]), str(dummy_expno)))

isDir = os.path.isdir(path_dummy_expno)
if isDir == False:
    os.mkdir(path_dummy_expno)
 #   os.mkdir(os.path.normpath(os.path.join(path_dummy_expno,"pdata")))

#######################################################
# Copy acquitision files from the InPhase dataset
if os.path.isfile(os.path.normpath(os.path.join(path_dummy_expno,'acqu'))) == False:
    for files in ['acqu','acqu2','acqus','acqu2s']:
        copyfile(os.path.normpath(os.path.join(current_dataset[3], current_dataset[0], str(expno_InPhase),files)),os.path.normpath(os.path.join(path_dummy_expno,files)))

#######################################################
# Run Processing with 1k in the direct dimension to perform the peak picking
RE([str(current_dataset[0]), str(expno_AntiPhase), "1", str(current_dataset[3])+"/"])
WR([str(current_dataset[0]), str(expno_AntiPhase), "1001", str(current_dataset[3])],"y")
RE([str(current_dataset[0]), str(expno_AntiPhase), "1001", str(current_dataset[3])+"/"])
PUTPAR('2 TDeff', str(1024))
XFB()

#######################################################
# Copy datasets into the dummy expno
for k, (ExpName, dest_procno, src_procno) in enumerate([
    (expno_InPhase, dummy_procno,   1), #InPhase
    (expno_AntiPhase, dummy_procno+1, 1), #AntiPhase
    (expno_AntiPhase, dummy_procno+2, 1001), #Antiphase 1k
    (expno_AntiPhase, dummy_procno+3, 1001), #Antiphase 1k
    (expno_InPhase, dummy_procno+4, 1), #Shifted
    ]):

    src_path = os.path.normpath(os.path.join(current_dataset[3],current_dataset[0],str(ExpName),"pdata",str(src_procno)))
    dest_path = os.path.normpath(os.path.join(path_dummy_expno,'pdata',str(dest_procno)))

    RE([str(current_dataset[0]), str(ExpName), str(src_procno), str(current_dataset[3])+"/"],show="y")

    # check spectrum dimension
    if GETPROCDIM() != 2:
        ERRMSG(message = "The spectrum to work with must have 2 dimensions.", title="Error", details=None, modal=1)
        EXIT()

    # copy_all_pdata_files(src_path,dest_path)
    if os.path.isdir(dest_path):
        rmtree(dest_path)

    os.makedirs(dest_path)
    copy_all_pdata_files(src_path,dest_path)

    # Create text file for each pdata
    file = os.path.normpath(os.path.join(dest_path,"title"))
    txt = ""
    txt += str(day)+'----'+str(current_time)+"\n"
    txt += "Analysis of Multiplet \n"
    txt += "\n"+Comment_Title
    if dest_procno == dummy_procno:
        txt += "\nExtracted from InPhase dataset : "+str(expno_InPhase) 
    if dest_procno == dummy_procno+1:
        txt += "\nExtracted from AntiPhase dataset : "+str(expno_AntiPhase) 
    if dest_procno == dummy_procno+2:
        txt += "\nExtracted from AntiPhase dataset : "+str(expno_AntiPhase) 
        txt += "\nProcessed with 1k points in the direct dimension for peak picking" 
    if dest_procno == dummy_procno+3:
        txt += "\nExtracted from AntiPhase dataset : "+str(expno_AntiPhase) 
        txt += "\nProcessed with 1k points in the direct dimension for peak picking" 
    if dest_procno == dummy_procno+4:
        txt += "\nShifted Spectrum" 
    if Wdw_ppm != [0,0,0,0]:
        txt += "\nRegion : "+str(Wdw_ppm[0])+ " ppm to "+str(Wdw_ppm[1])+ " ppm // "+str(Wdw_ppm[2])+ " ppm to "+str(Wdw_ppm[3])+ " ppm "  

    f = open(file, "w")
    f.write(txt)
    f.close()

Path                = str(current_dataset[3])     #1
Data_Folder         = str(inputs[0])              #2
InPhase_DataSet     = expno_InPhase               #3
AntiPhase_DataSet   = expno_AntiPhase             #4  
Dummy_DataSet       = dummy_expno                 #5
Dummy_procno        = dummy_procno                #6

print("Start Analysis Script")
cmd = tr(os.path.normpath(os.path(options["python_path"])))+" "+str(Analysis_Script_Path)+" %s %s %s %s %s %s %s %s %s %s %s %s %s"  % (Path, Data_Folder, InPhase_DataSet, AntiPhase_DataSet, Dummy_DataSet, Dummy_procno, Wdw_ppm[0], Wdw_ppm[1], Wdw_ppm[2], Wdw_ppm[3], pp_th, Box_pts[0], Box_pts[1])
subprocess.Popen(cmd, shell=True).wait()
print('End Analysis Script')
EXIT()

