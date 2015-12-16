#!/usr/bin/python

# Command line LCMODEL analysis of Siemens data
#    File Name:  process_fmri.py
#
#    NOTES - 2014/02/27 - Wayne Lee - Initial creation build off rdatools.py created by Dallas Card
#
#  Crontab for lcmodel on sagan (put on a single line)
#    30 23 * * * for dir_mrs in /mnt/ResearchMR/GE_Research/mrs_xfr/MRS_ONLY/*; do 
#           /usr/local/bin/LCM_tools/LCM_siemens.py -d -v -m --logfile=LogFile.txt --dir_backup /mnt/ResearchMR/ResearchPACS ${dir_mrs}; done


from optparse import OptionParser
import subprocess
import os.path
import numpy
import struct
import sys
import datetime
import string

program_name = 'LCM_siemens.py'
valid_chars = '-_%s%s' % (string.ascii_letters, string.digits)
lut_fieldstrength = {};
lut_fieldstrength['MRC35421'] = '3t'
lut_fieldstrength['MRC26247'] = '1.5t'

lut_scanner_type = {};
lut_scanner_type['MRC35421'] = 'siemens'
lut_scanner_type['MRC26247'] = 'siemens'


def run_cmd(sys_cmd, options):
# one line call to output system command and control debug state
    if options.verbose:
        print " > " + sys_cmd
    if not options.debug:
        p = subprocess.Popen(sys_cmd, stdout = subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        output, errors = p.communicate()
        return output, errors
    else:
        return '',''        

# Check if file exists before running
def check_and_run(sys_cmd, dir_out, prefix_out, file_type, options):
    # if file exists
    if os.path.exists( '%s/%s%s' % (dir_out, prefix_out, file_type)):
        # Check if it's clobbering time
        if options.clobber:    
            print " > Overwriting file - %s/%s%s" % (dir_out, prefix_out, file_type)
            clobber_cmd = 'rm %s/%s%s' % (dir_out, prefix_out, file_type)
            output, errors = run_cmd(clobber_cmd, options)
            
            output, errors = run_cmd(sys_cmd, options)
            return output, errors
        else:
            print ' > %s > Output file already exists ' % (sys_cmd,)
            return '',''
    else:     # file doesn't exist
        output, errors = run_cmd(sys_cmd, options)
        # Check after program run for output file's existance if not in debug mode
        if (not options.debug) and (not os.path.exists( '%s/%s%s' % (dir_out, prefix_out, file_type))):
            print '%%% ERROR - Function call failed - %s' % (sys_cmd,) 
            print ' > Function call output:'
            print output, errors
            # raise SystemExit, 'ERROR - Function call failed - %s' % (sys_cmd,) 
        return output, errors

        

# FUNCTION TO CHECK AND MAKE A DIRECTORY IF NECESSARY
def check_dir(dir, options):
    if not os.path.exists(dir):
        run_cmd('mkdir ' + dir, options)
    else:
        print ' > mkdir ' + dir + ' > Directory already exists'
    return dir
    
# import_rda(filename)
# This function reads in a .rda file containing spectroscopy data from the Siemens scanner.
# It parses the header, places the relevant information into a dictionary, and
# returns the header and the data (fid) as an array of complex values
def import_rda(filename):
    # define empty arrays
    header_names = []
    header_vals = []
    
    # open the file
    file = open(filename)

    # read header and first line
    line = file.readline()
    line = file.readline()
    while (not ("End of header" in line)):
        parts = line.split(':')
        header_names.append(parts[0])
        header_vals.append(parts[1].rstrip('\r\n'))
        line = file.readline()

    # make a dictionary of the header values
    header = dict(zip(header_names, header_vals))

    # get the number of samples
    num_samples = int(header.get("VectorSize"))

    # make a variable to hold the fid, and read it in from the file
    fid = numpy.zeros(num_samples, dtype=numpy.complex)
    for i in range(num_samples):
        real_part = struct.unpack('d', file.read(8))
        imag_part = struct.unpack('d', file.read(8))
        fid[i] = real_part[0] - imag_part[0] * 1j

    # close the file
    file.close()
    
    # return the fid and dictionary of header information
    return header, fid
    

# get_params(header)
# Pull and format relevant parameters from the header
def get_params(header):
    params = {};
    
    list_string = ('SeriesDate','PatientName','SeriesNumber','SeriesDescription','StationName', \
        'SeriesTime', 'PatientID')
    for name_string in list_string:
        params[name_string] = ''.join(c for c in header.get(name_string) if c in valid_chars)
 
    list_int = ('TE','TR','VectorSize','InstanceNumber')
    for name_int in list_int:
        params[name_int] = int(float(header.get(name_int).strip(' ')))
        
    list_float = ('DwellTime','MRFrequency')
    for name_float in list_float:
        params[name_float] = float(header.get(name_float).strip(' ')) 

    params['SeriesTime'] = params['SeriesTime'][0:2] + ':' + params['SeriesTime'][2:4]
    params['FOV_read'] = float(header.get('VOIReadoutFOV').strip(' '))        # mm
    params['FOV_phase'] = float(header.get('VOIPhaseFOV').strip(' '))         # mm
    params['FOV_slice'] = float(header.get('VOIThickness').strip(' '))        # mm
    params['Voxel_Volume'] = params['FOV_read'] * params['FOV_phase'] * params['FOV_slice'] / 1000   # ml
    params['DwellTime'] = params['DwellTime'] / 1000000 # s
    params['BW'] = 1 / params['DwellTime']              # Hz
    params['Hz_per_ppm'] = params['MRFrequency']
    
    return params

    
# check _params(hdr_met, hdr_h2o)
# Checks metabolite and hdr parameters to ensure that they match and "should" be analyzed together
def check_params(hdr_met, hdr_h2o):
    list_params_to_check = ('SeriesDate','PatientName','StationName','PatientID','TE','TR','DwellTime', \
        'FOV_read','FOV_phase','FOV_slice')
    params_pass = 1
    for curr_param in list_params_to_check:
        if hdr_met[curr_param] != hdr_h2o[curr_param]:
            params_pass = 0
            print "*** ERROR - met and h2o parameters do not match [%s: %s != %s] *** " % \
                (curr_param, hdr_met[curr_param], hdr_h2o[curr_param])
            
    if not params_pass:
        sys.exit()
    
    
# get_basis(params)
# Determine the basis set to use
def get_basis(params, options):

    # SIEMENS BASIS SETS
    if params['TE'] < 70.5:
        basis_set = 'gamma_press_te30_64mhz_083.basis' 
    else:
        basis_set = 'gamma_press_te144_64mhz_083.basis' 
    file_basis = '%s/basis-sets/siemens/%s' % \
        (options.dir_lcmodel, basis_set)

    return file_basis;
    

# create_control(dir_series, params, options)
# create control file for analysis
def create_control(dir_series, title_LCM, params, options):
    if not options.debug:
        f_control = open( dir_series + '/control' ,'w+')
        f_control.write(' $LCMODL\n')
        f_control.write(" title= '%s'\n" % (title_LCM,))

        f_control.write(' sddegz= 3.\n')  # Std dev for expected 0-order phase corr
        f_control.write(' sddegp= 1.\n')  # Std dev for expected 1st-order phase corr
        f_control.write(" savdir= '%s'\n" % (dir_series,))
        f_control.write(' ppmst= 4.0\n')   # start of spectra 
        f_control.write(' ppmend= 0.2\n')  # end of spectra 
        f_control.write(' nunfil= %d\n' % (params['VectorSize'],))
        f_control.write(' ltable= 7\n')     # 7= create a table output
        f_control.write(' lps= 8\n')        # 8 = create ps output
        f_control.write(' lcsv= 11\n')      # 11 = create csv output
        f_control.write(' hzpppm= %.4e\n' % (params['Hz_per_ppm'],) )
        f_control.write(" filtab= '%s/table'\n" % (dir_series,))
        f_control.write(" filraw= '%s/met/RAW'\n" % (dir_series,))
        f_control.write(" filps= '%s/ps'\n" % (dir_series,))
        f_control.write(" filh2o= '%s/h2o/RAW'\n" % (dir_series,))
        f_control.write(" filcsv= '%s/spreadsheet.csv'\n" % (dir_series,))
        f_control.write(" filbas= '%s'\n" % (get_basis(params,options),))
        f_control.write(' echot= %d\n' % (int(params['TE']),))
        f_control.write(' dows= T\n')       # Do water scaling
        f_control.write(' doecc= T\n')      # Do Eddy Current Correction
        f_control.write(' deltat= %.4e\n' % (params['DwellTime'],))
        f_control.write(' $END\n')
        f_control.close()

    
if __name__ == '__main__' :   
    usage = "usage: [DIRECTORY MODE] "+program_name+" <options> dir_rda \n" + \
            "                 dir_rda must contain a single pair of rda files \n" + \
            "                     XXX_met.rda - Metabolite acquisition (suffix must be exact)\n" + \
            "                     XXX_h2o.rda - Water reference acquisitions (suffix must be exact) \n" + \
            "   or  [FILE MODE] "+program_name+" <options> fname_rda_met fname_rda_h2o \n" + \
            "   or  "+program_name+" -h";
            
    parser = OptionParser(usage=usage)
    parser.add_option("-c","--clobber", action="store_true", dest="clobber",
                        default=0, help="overwrite output file")
    parser.add_option("-d","--debug", action="store_true", dest="debug",
                        default=0, help="Run in debug mode")
    parser.add_option("-v","--verbose", action="store_true", dest="verbose",
                        default=0, help="Verbose output")
    parser.add_option("-m","--move_original", action="store_true", dest="move_original",
                        default=0, help="Move original files")
    parser.add_option("--dir_backup", type='string', dest="dir_backup",
                        default='', help="If a directory is supplied, base directory for backing up rda file [%default]")

    # read in arguments
    (options, args) = parser.parse_args()
    
    # Validate input arguments
    if len(args) == 0:
        print usage
        sys.exit()
    elif len(args) == 1:
        fname_rda = args[0]
    else:
        print "\n*** ERROR - Incorrect number of arguments! (1 expected, %s received) ***\n" % (len(args),)
        print usage
        sys.exit()
        
    # Load rda files and split into headers and data
    hdr_rda, fid_rda = import_rda(fname_rda)   
    # print hdr_rda
    # Extract and validate relevant parameters from headers
    params = get_params(hdr_rda)

    # Backup RDA 
    check_dir(options.dir_backup,options)
    options.dir_backup = '%s/%s-%s' % (options.dir_backup, params['SeriesDate'], params['PatientName'])
    check_dir(options.dir_backup,options)
    options.dir_backup = '%s/mrs' % (options.dir_backup,)
    check_dir(options.dir_backup,options)
    
    if options.move_original:
        cmd_handle_original = 'mv -u'
    else:
        cmd_handle_original = 'cp -au' 
    
    fname_output = '%s-%s-%d' % ( params['SeriesNumber'], params['SeriesDescription'], params['InstanceNumber'])
    cmd_backup = '%s %s %s/%s.rda' % (cmd_handle_original, fname_rda , options.dir_backup, fname_output)
    check_and_run(cmd_backup, options.dir_backup, fname_output, '.rda', options)
   

   
