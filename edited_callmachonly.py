import re, os, glob, subprocess, shutil

 
# get current working dir
WORK_DIR = os.getcwd()
#########################################################################
# PARAMETERS, I.E. INPUTS TO KNIFE CALL; need to add these here
#########################################################################

# should you use the toy indel indices; use them for testing; generally this should be False
use_toy_indel = True

# dataset_name CANNOT HAVE ANY SPACES IN IT
import argparse
parser = argparse.ArgumentParser()
parser.add_argument("--dataset", help="name of dataset-NO SPACES- will be used for naming files")
# Should fix this later: not so good, no time now; get rid of dataset argument
# and search for dataset name
parser.add_argument("--knifeoutputtarball", help="path to tarred and zipped file of knife output directory, the tarball should unpack with the name of the directory the same as the dataset name")
parser.add_argument("--runid", help="name of dataset-NO SPACES- will be used for naming files")


args = parser.parse_args()
if args.dataset:
    dataset_name = args.dataset
else:
    dataset_name = "nodatasetname"

knife_output_tarball = args.knifeoutputtarball

if args.runid:
    run_id = args.runid
else:
    run_id = "norunid"
    

    
# e.g. dataset_name = "4ep18" 

report_directory_name = "circReads"


# mode = "skipDenovo"
# read_id_style= "appended"
# junction_overlap =  8
# ntrim= 40

# Not really used, just doing so it mimics test Data call
# logstdout_from_knife = "logofstdoutfromknife"


logfile = WORK_DIR + "/logmachonly" + dataset_name + run_id +".txt"

with open(logfile, 'w') as ff:
    ff.write(WORK_DIR)
    ff.write('\n\n\n')

    
# main directory to be used when running the knife:
#knifedir = "/srv/software/knife/circularRNApipeline_Standalone"
cirpipedir = "/srv/software/knife/circularRNApipeline_Standalone"

for file in os.listdir(cirpipedir):
	src_file = os.path.join(cirpipedir, file)
	dst_file = os.path.join(WORK_DIR, file)
	shutil.move(src_file, dst_file) 

knifedir = WORK_DIR

#########################################################################
#
# tar and unpack the knife output directory
#   check that there is then a directory called dataset_name
#   this is a hack for now
#
#########################################################################

os.chdir(WORK_DIR)
with open(logfile, 'a') as ff:
    ff.write("\nAbout to try to unpack the tarfile " + knife_output_tarball + "\n")
try:
    fullcall = "tar -xvzf " + knife_output_tarball
    with open(logfile, 'a') as ff:
        subprocess.check_call(fullcall, stderr=ff, stdout = ff, shell=True)
except:
    with open(logfile, 'a') as ff:
        ff.write("\nError in unpacking the tarfile " + knife_output_tarball + "\n")

datadir = os.path.join(WORK_DIR,dataset_name)
if not os.path.isdir(datadir):
    with open(logfile, 'a') as ff:
        ff.write('ERROR: no directory\n' + datadir + '\nMake sure input of dataset matches dataset folder name unpacked from tarball.\n All later work is suspect.\n')
        
        

#########################################################################
#
# Note: should have unpacked tar of knife output results before doing
#   next part.
#
#########################################################################


###Run MACHETE
#Nathaniel Watson
#05-26-2016 

CIRCPIPE_DIR = os.path.join(WORK_DIR,dataset_name)
if not os.path.isdir(CIRCPIPE_DIR):
    with open(logfile, 'a') as ff:
        ff.write('Problem: no directory\n' + CIRCPIPE_DIR + '\nMaking one.')
    os.mkdir(CIRCPIPE_DIR)
CIRCREF = os.path.join(knifedir,"index")
# CIRCREF = os.path.join(WORK_DIR,dataset_name,report_directory_name,"index")
if not os.path.isdir(CIRCREF):
    with open(logfile, 'a') as ff:
        ff.write('No directory\n' + CIRCREF + '\nNot so surprising.\nMaking one.')
    os.mkdir(CIRCREF)
MACH_OUTPUT_DIR = os.path.join(WORK_DIR,"mach")
os.mkdir(MACH_OUTPUT_DIR)
EXONS = os.path.join(WORK_DIR,"HG19exons")

if use_toy_indel:
    REG_INDEL_INDICES = os.path.join(WORK_DIR,"toyIndelIndices") #test indices for faster runs
else:
    REG_INDEL_INDICES = os.path.join(WORK_DIR,"IndelIndices")

#########################################################################
#
# unpack HG19exons.tar.gz and IndelIndices.tar.gz
#
#########################################################################

os.chdir(WORK_DIR)

# Note that you have to choose a directory to start in, to unpack there
def unpack_tarball_with_checking(tarball, thislogfile, thisworkdir):
    os.chdir(thisworkdir)
    if os.path.isfile(tarball):
        with open(thislogfile, 'a') as ff:
            ff.write("\nAbout to try to unpack the tarfile " + tarball + "\n")
        try:
            fullcall = "tar -xvzf " + tarball
            with open(thislogfile, 'a') as ff:
                subprocess.check_call(fullcall, stderr=ff, stdout = ff, shell=True)
        except:
            with open(thislogfile, 'a') as ff:
                ff.write("\nError in unpacking the tarfile " + tarball + " \n")
    else:
        with open(thislogfile, 'a') as ff:
            ff.write("\nERROR: No tarball found called " + tarball + " \n")

unpack_tarball_with_checking(tarball = "HG19exons.tar.gz", thislogfile=logfile, thisworkdir=WORK_DIR)
if use_toy_indel:
    unpack_tarball_with_checking(tarball = "toyIndelIndices.tar.gz", thislogfile=logfile, thisworkdir=WORK_DIR)
else:
    unpack_tarball_with_checking(tarball = "IndelIndices.tar.gz", thislogfile=logfile, thisworkdir=WORK_DIR)
            
            

#########################################################################
#
# move reference libraries output by KNIFE to directory CIRCREF that
#   contains hg19_genome, hg19_transcriptome, hg19_junctions_reg and
#   hg19_junctions_scrambled bowtie indices
#
#########################################################################
    
 
# first have to 
#   change file names and mv them to the right directories


# Input file names are in an unusual format so they are easy to select when doing a run on
#   seven bridges. They should start in the home directory, as copies, because
#   they are entered as stage inputs. They start with the prefix "infile"
#   Move them to the directory where MACHETE
#   expects them to be, then change their names.


# knifedir = "/srv/software/knife/circularRNApipeline_Standalone"
knifedir = WORK_DIR

    
    
targetdir_list = [knifedir + "/index", knifedir + "/denovo_scripts", knifedir + "/denovo_scripts/index"]
    
# check that all the subdirectories are there, as they should be!
# these should have been made beforehand, and should have appropriate files in them:
# files starting with infilebt2 and infilefastas should be in /index
#   and the prefixes "infilebt2" and "infilefastas" should be removed
# files starting with infilegtf should be in /denovo_scripts
#   and the prefix "infilegtf" should be removed
# files starting with infilebt1 should be in /denovo_scripts/index
#   and the prefix "infilebt1" should be removed

thisdir = targetdir_list[0]
if not os.path.exists(thisdir):
    raise ValueError("Error: directory " + thisdir + " does not exist.")

thisdir = targetdir_list[1]
if not os.path.exists(thisdir):
    raise ValueError("Error: directory " + thisdir + " does not exist.")

thisdir = targetdir_list[2]
if not os.path.exists(thisdir):
    raise ValueError("Error: directory " + thisdir + " does not exist.")

# cd into the knife directory; should not really be necessary
os.chdir(knifedir)

with open(logfile, 'a') as ff:
    ff.write('\n\n\n')
    ff.write('Listing files in knifedir ' + knifedir) 
    
with open(logfile, 'a') as ff:
    ff.write('\n\n\n')
    subprocess.check_call(["ls", "-R"], stdout=ff)
    ff.write('\n\n\n')


# Did the below in past, but now just change name of CIRCREF in run.py
    
# # Should still be in working dir now, but in case not
# os.chdir(WORK_DIR)

# # Note that files have prefixes like infilebt1 and infilebt2 b/c the extra
# #  bt1 and bt2 told the knife program where to send them

# prefix = "infile"
# globpattern = prefix + "*"
# matching_files = glob.glob(globpattern)
# if (len(matching_files)>= 1):
#     for thisfile in matching_files:
#         fullpatholdfile = WORK_DIR + "/" + thisfile
#         fullpathnewfile = CIRCREF + "/" + re.sub(pattern=prefix, repl="", string= thisfile)
#         with open(logfile, 'a') as ff:
#             ff.write('About to do mv '+ fullpatholdfile + ' ' + fullpathnewfile + '\n')
#         with open(logfile, 'a') as ff:
#             subprocess.check_call(["mv", fullpatholdfile, fullpathnewfile], stderr=ff, stdout=ff)
        

# # run test of knife
# # sh completeRun.sh READ_DIRECTORY complete OUTPUT_DIRECTORY testData 8 phred64 circReads 40 2>&1 | tee out.log

# try:
#     with open(logfile, 'a') as ff:
#         ff.write('\n\n\n')
#         # changing so as to remove calls to perl:
#         subprocess.check_call("sh completeRun.sh " + WORK_DIR + " " + read_id_style + " " + WORK_DIR + " " + dataset_name + " " + str(junction_overlap) + " " + mode + " " + report_directory_name + " " + str(ntrim) + " 2>&1 | tee " + logstdout_from_knife , stdout = ff, shell=True)
#         # original test call:
#         # subprocess.check_call("sh completeRun.sh " + WORK_DIR + " complete " + WORK_DIR + " testData 8 phred64 circReads 40 2>&1 | tee outknifelog.txt", stdout = ff, shell=True)
# except:
#     with open(logfile, 'a') as ff:
#         ff.write('Error in running completeRun.sh')



# datadirlocation = WORK_DIR + "/" + dataset_name  

#############################################################################
#
# Now run the machete
#
#############################################################################




MACH_DIR = "/srv/software/machete"
MACH_RUN_SCRIPT = os.path.join(MACH_DIR,"run.py")

cmd = "python {MACH_RUN_SCRIPT} --circpipe-dir {CIRCPIPE_DIR} --output-dir {MACH_OUTPUT_DIR} --hg19Exons {EXONS} --reg-indel-indices {REG_INDEL_INDICES} --circref-dir {CIRCREF}".format(MACH_RUN_SCRIPT=MACH_RUN_SCRIPT,CIRCPIPE_DIR=CIRCPIPE_DIR,MACH_OUTPUT_DIR=MACH_OUTPUT_DIR,EXONS=EXONS,REG_INDEL_INDICES=REG_INDEL_INDICES,CIRCREF=CIRCREF)

with open(logfile, 'a') as ff:
        ff.write('\n\n\nAbout to run run.py\n')
        ff.write('\n\n\n')

with open(logfile, 'a') as ff:
        retcode = subprocess.call(cmd,shell=True, stdout=ff, stderr=ff)
#        popen = subprocess.check_call(cmd,shell=True, stdout=ff, stderr=ff)

## Write either YES*.error.in.subprocess.txt with 1 or
##     NO*.error.in.subprocess.txt with 0
## This allows you to see quickly in the results on Seven Bridges if there is 
##     an error in the subprocess call in run.py
## File ends in error.in.subprocess.txt; the name changes depending on the value.

errorfiles = [os.path.join(WORK_DIR, x + ".error." + dataset_name + run_id + '.is.error.in.subprocess.txt')  for x in ['YES','NO']]

if retcode:
    with open(errorfiles[0], 'w') as ff:
        ff.write('1')
else:
    with open(errorfiles[1], 'w') as ff:
        ff.write('0')


        

os.chdir(WORK_DIR)


with open(logfile, 'a') as ff:
    ff.write('\n\n\nListing working directory but not recursively.\n\n\n')
    
with open(logfile, 'a') as ff:
    subprocess.check_call("ls -alh", stdout=ff, stderr=ff, shell=True)
    
os.chdir(MACH_OUTPUT_DIR)


with open(logfile, 'a') as ff:
    ff.write('\n\n\nListing machete ouput directory recursively.\n\n\n')
    
with open(logfile, 'a') as ff:
    subprocess.check_call("ls -alRh", stdout=ff, stderr=ff, shell=True)
    
with open(logfile, 'a') as ff:
    ff.write('\n\n\nMasterError.txt should start here.\n\n\n')

os.chdir(WORK_DIR)
    
fullpatholderrorfile = MACH_OUTPUT_DIR + "/MasterError.txt"
if os.path.isfile(fullpatholderrorfile):
    subprocess.check_call("cat " + fullpatholderrorfile + " >> " + logfile, shell=True)
else:
    with open(logfile, 'a') as ff:
        ff.write('\n\n\nNo MasterError.txt file found.\n\n\n')

# tar everything in mach_output_dir/reports 
#   Tar if it exists
os.chdir(WORK_DIR)
mach_output_reports_dir = os.path.join(MACH_OUTPUT_DIR,"reports")
if os.path.isdir(mach_output_reports_dir):
    try:
        fullcall = "tar -cvzf " + dataset_name + run_id + "machreportsout.tar.gz -C " + MACH_OUTPUT_DIR + " reports"  
        with open(logfile, 'a') as ff:
            subprocess.check_call(fullcall, stderr=ff, stdout = ff, shell=True)
    except:
        with open(logfile, 'a') as ff:
            ff.write("\nError in tarring the machete output report files in the " + mach_output_reports_dir + " directory\n")
else:
    with open(logfile, 'a') as ff:
        ff.write("\nNo directory of machete output reports called " + mach_output_reports_dir + ", but expected that there was one.\n")
    

# tar everything in mach_output_dir/err_and_out if the directory exists
os.chdir(WORK_DIR)
mach_err_dir = os.path.join(MACH_OUTPUT_DIR,"err_and_out")
if os.path.isdir(mach_err_dir):
    try:
        fullcall = "tar -cvzf " + dataset_name + run_id + "macherrout.tar.gz -C " + MACH_OUTPUT_DIR + " err_and_out"  
        with open(logfile, 'a') as ff:
            subprocess.check_call(fullcall, stderr=ff, stdout = ff, shell=True)
    except:
        with open(logfile, 'a') as ff:
            ff.write("\nError in tarring the machete err_and_out directory, namely the " + mach_err_dir + " directory\n")
else:
    with open(logfile, 'a') as ff:
        ff.write("\nNo directory of machete errors and output called " + mach_err_dir + ", but expected that there was one.\n")

# tar everything in mach_output_dir/BadFJ_ver2 if the directory exists- changing name on jun 20
os.chdir(WORK_DIR)
mach_err_dir = os.path.join(MACH_OUTPUT_DIR,"BadFJ_ver2")
if os.path.isdir(mach_err_dir):
    try:
        fullcall = "tar -cvzf " + dataset_name + run_id + "machbadfjver2out.tar.gz -C " + MACH_OUTPUT_DIR + " BadFJ_ver2"  
        with open(logfile, 'a') as ff:
            subprocess.check_call(fullcall, stderr=ff, stdout = ff, shell=True)
    except:
        with open(logfile, 'a') as ff:
            ff.write("\nError in tarring the machete BadFJ_ver2 directory, namely the " + mach_err_dir + " directory\n")
else:
    with open(logfile, 'a') as ff:
        ff.write("\nNo directory of machete errors and output called " + mach_err_dir + ", but expected that there was one.\n")


# tar everything in mach_output_dir/x if the directory exists

def tar_subdirectory_of_mach_output_dir(thisdir, text_for_naming, MACH_OUTPUT_DIR, WORK_DIR, dataset_name, run_id, logfile):
    os.chdir(WORK_DIR)
    thisfulldir = os.path.join(MACH_OUTPUT_DIR,thisdir)
    if os.path.isdir(thisfulldir):
        try:
            fullcall = "tar -cvzf " + dataset_name + run_id + "mach" + text_for_naming + "out.tar.gz -C " + MACH_OUTPUT_DIR + " " + thisdir  
            with open(logfile, 'a') as ff:
                subprocess.check_call(fullcall, stderr=ff, stdout = ff, shell=True)
        except:
            with open(logfile, 'a') as ff:
                ff.write("\nError in tarring the machete " + thisdir + " directory, namely the " + thisfulldir + " directory\n")
    else:
        with open(logfile, 'a') as ff:
            ff.write("\nNo directory within machete output directory called " + thisfulldir + ", but expected that there was one.\n")

tar_subdirectory_of_mach_output_dir(thisdir = "GLM_classInput", text_for_naming = "glmclassinput", MACH_OUTPUT_DIR=MACH_OUTPUT_DIR, WORK_DIR=WORK_DIR, dataset_name=dataset_name, run_id=run_id, logfile=logfile)

tar_subdirectory_of_mach_output_dir(thisdir = "reports", text_for_naming = "reports", MACH_OUTPUT_DIR=MACH_OUTPUT_DIR, WORK_DIR=WORK_DIR, dataset_name=dataset_name, run_id=run_id, logfile=logfile)

tar_subdirectory_of_mach_output_dir(thisdir = "BadFJ", text_for_naming = "badfj", MACH_OUTPUT_DIR=MACH_OUTPUT_DIR, WORK_DIR=WORK_DIR, dataset_name=dataset_name, run_id=run_id, logfile=logfile)

# Now also get tar file of glm reports from knife for convenient reading together with
#  machete output.

os.chdir(WORK_DIR)

glmdirlocation = os.path.join(WORK_DIR, dataset_name, "circReads")

knifeglmreportstarfile = dataset_name + run_id + "knifeglmreportfilesout.tar.gz"

fullcall = "tar -cvzf " + knifeglmreportstarfile + " -C " + glmdirlocation + " glmReports"

with open(logfile, 'a') as ff:
    subprocess.check_call(fullcall, stdout=ff, stderr=ff, shell=True)

# Now copy StemList.txt and MasterError.txt to
#  files with better names so that they have unique names for each run
stemoutfile = WORK_DIR + "/StemList" + dataset_name + run_id +".txt"
stemorigfile = os.path.join(MACH_OUTPUT_DIR, "StemList.txt")
fullcall = "cp " + stemorigfile + " " + stemoutfile 
with open(logfile, 'a') as ff:
    subprocess.check_call(fullcall, stdout=ff, stderr=ff, shell=True)

mastererroroutfile = WORK_DIR + "/MasterError" + dataset_name + run_id +".txt"
mastererrororigfile = os.path.join(MACH_OUTPUT_DIR, "MasterError.txt")
fullcall = "cp " + mastererrororigfile + " " + mastererroroutfile
with open(logfile, 'a') as ff:
    subprocess.check_call(fullcall, stdout=ff, stderr=ff, shell=True)
