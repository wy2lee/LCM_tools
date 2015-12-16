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
 
    list_int = ('TE','TR','VectorSize')
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
    if lut_fieldstrength[ params['StationName']] == '3t':
        if params['TE'] < 27.5:
            basis_set = 'press_te25_3t_gsh_v3.basis'
        elif params['TE'] < 32.5:
            basis_set = 'press_te30_3t_gsh_v3.basis'
        elif params['TE'] < 70:
            basis_set = 'gamma_press_te35_128mhz_087.basis'
        elif params['TE'] < 138:
            basis_set = 'press_te135_3t_gsh_v3.basis'
        else:
            basis_set = 'gamma_press_te144_128mhz_087.basis'
    elif lut_fieldstrength[ params['StationName']] == '1.5t':
        if params['TE'] < 32.5:
            basis_set = 'press_te30_64mhz_gsh_v3.basis'
        elif params['TE'] < 70:
            basis_set = 'press_te35_64mhz_gsh_v3.basis'
        elif params['TE'] < 160:
            basis_set = 'press_te144_64mhz_gsh_v3.basis'
        else:
            basis_set = 'press_te288_64mhz_gsh_v3.basis'

    print lut_fieldstrength[params['StationName']]
    file_basis = '%s/basis-sets/%s/%s' % \
        (options.dir_lcmodel, lut_fieldstrength[params['StationName']], basis_set)

    return file_basis;
    

# create_control(dir_series, params, options)
# create control file for analysis
def create_control(dir_series, title_LCM, params, options):
    if not options.debug:
        f_control = open( dir_series + '/control' ,'w+')
        f_control.write(' $LCMODL\n')
        f_control.write(" title= '%s'\n" % (title_LCM,))

#        f_control.write(' sddegz= 3.\n')  # Std dev for expected 0-order phase corr
#        f_control.write(' sddegp= 1.\n')  # Std dev for expected 1st-order phase corr
        f_control.write(" savdir= '%s'\n" % (dir_series,))
        f_control.write(' ppmst= 4.0\n')   # start of spectra 
        f_control.write(' ppmend= 1.0\n')  # end of spectra 
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
    parser.add_option("--no_h2o", action="store_true", dest="no_h2o",
                        default=0, help="No h2o scan")
    parser.add_option("--printer", type='string', dest="printer",
                        default='', help="If supplied, print ps output to this printer")
    parser.add_option("--dir_backup", type='string', dest="dir_backup",
                        default='', help="If a directory is supplied, base directory for backing up rda file")
    parser.add_option("--dir_lcmodel", type='string', dest="dir_lcmodel",
                        default='/home/lcmodel/.lcmodel', help="Directory of lcmodel install [%default]")
    parser.add_option("--dir_output", type='string', dest="dir_output",
                        default='/data/LCM_DATA', help="Base output directory [%default]")
    parser.add_option("--logfile", type='string', dest="logfile",
                        default='', help="If supplied, update logfile in <dir_output>/<scanner>/")

    # read in arguments
    (options, args) = parser.parse_args()
    
    # Validate input arguments
    if len(args) == 0:
        print usage
        sys.exit()
    elif (len(args) == 1) and (options.no_h2o==0):
        dir_input_rda = args[0]
	if not os.path.exists(dir_input_rda):
            print "\n*** No directories in path *** \n"
            sys.exit()
        ldir_input_rda = os.listdir(dir_input_rda)
        count_met = 0
        count_h2o = 0
        # check for single h2o and met rda files
        for fname_temp in ldir_input_rda:
            if fname_temp.find('met.rda') > -1:
                fname_rda_met = '%s/%s' % (dir_input_rda,fname_temp)
                count_met = count_met + 1
            if fname_temp.find('h2o.rda') > -1:
                fname_rda_h2o = '%s/%s' % (dir_input_rda,fname_temp)
                count_h2o = count_h2o + 1
        if count_met != 1:
            print "\n*** ERROR - Too many/few XXX_met.rda files, there can only be one! [%d found] *** \n" % \
                (count_met,)
        if count_h2o != 1:
            print "\n*** ERROR - Too many/few XXX_h2o.rda files, there can only be one! [%d found] *** \n" % \
                (count_h2o,)
        if (count_met != 1) or (count_h2o != 1):
            sys.exit()
    elif (len(args) == 1) and (options.no_h2o ==1):  # No water reference
        if os.path.isdir(args[0]):
            dir_input_rda = args[0]
            ldir_input_rda = os.listdir(dir_input_rda)
            count_met = 0
            for fname_temp in ldir_input_rda:
                if fname_temp.find('met.rda') > -1:
                    fname_rda_met = '%s/%s' % (dir_input_rda,fname_temp)
                    count_met = count_met + 1
            if count_met != 1:
                print "\n*** ERROR - Too many/few XXX_met.rda files, there can only be one! [%d found] *** \n" % \
                    (count_met,)
                sys.exit()
            fname_rda_h2o = "NA"
        else:
            fname_rda_met = args[0]
            fname_rda_h2o = "NA"
    elif len(args) == 2:
        fname_rda_met, fname_rda_h2o = args
    else:
        print "\n*** ERROR - Incorrect number of arguments! (1 or 2 expected, %s received) ***\n" % (len(args),)
        print usage
        sys.exit()
        
    # Load rda files and split into headers and data
    hdr_met, fid_met = import_rda(fname_rda_met)   
    if not options.no_h2o:
	hdr_h2o, fid_h2o = import_rda(fname_rda_h2o)

    # Extract and validate relevant parameters from headers
    params = get_params(hdr_met)
    if not options.no_h2o:
	params_h2o = get_params(hdr_h2o)
	check_params(params, params_h2o)

    # Create LCM title for this analysis
    title_LCM = '%s - %s (%s) %s-%s [TR/TE=%d/%d; %4.3fml]' % \
        (params['SeriesDate'], params['PatientName'],  params['PatientID'], \
        params['SeriesNumber'], params['SeriesDescription'], \
        params['TR'], params['TE'], params['Voxel_Volume'])
    if options.verbose:
        print '%%% PROCESSING - ' + title_LCM
    
    # Prepare output directories
    if options.verbose:
        print '%%% Prepping output directories'

    name_study = '%s-%s' % (params['SeriesDate'], params['PatientName'])
    name_series = '%s-%s' % (params['SeriesNumber'], params['SeriesDescription'])
    if not options.no_h2o:
        name_h2o = '%s-%s' % (params_h2o['SeriesNumber'], params['SeriesDescription'])

    dir_study = '%s/%s/processed/%s' % \
        (options.dir_output, params['StationName'], name_study)
    check_dir(dir_study, options)

    dir_series = '%s/%s' % (dir_study, name_series)
    check_dir(dir_series, options)

    dir_met = '%s/met' % (dir_series,)
    check_dir(dir_met, options)

    if not options.no_h2o:
        dir_h2o = '%s/h2o' % (dir_series,)
        check_dir(dir_h2o, options)
    
    # Backup RDA to alternate directory if necessary
    # Matches input format (ie. within a sub-directory or just as two files)
    if not options.dir_backup == '':
        if options.verbose:
            print '%%% Creating RDA backup'
        check_dir(options.dir_backup,options)
        
        if len(args)==1:
            options.dir_backup = '%s/%s' % (options.dir_backup, name_study)
            check_dir(options.dir_backup,options)
            options.dir_backup = '%s/mrs' % (options.dir_backup)
            check_dir(options.dir_backup,options)
            options.dir_backup = '%s/%s' % (options.dir_backup, name_series)
            check_dir(options.dir_backup,options)
        
        cmd_backup_met = 'cp -au %s %s/%s_met.rda' % (fname_rda_met , options.dir_backup, name_series)
        run_cmd(cmd_backup_met, options)
        if not options.no_h2o:
            cmd_backup_h2o = 'cp -au %s %s/%s_h2o.rda' % (fname_rda_h2o , options.dir_backup, name_h2o)
            run_cmd(cmd_backup_h2o, options)    

        
    # Make a copy of RDAs in LCM output directory
    if options.verbose:
        print '%%% Making a copy of RDA files in LCM output directory'
    if options.move_original:
        cmd_handle_original = 'mv -f'
    else:
        cmd_handle_original = 'cp -f' 
    
    cmd_cp_met = '%s %s %s/%s_met.rda' % (cmd_handle_original, fname_rda_met , dir_series, name_series)
    run_cmd(cmd_cp_met, options)

    if not options.no_h2o:   
        cmd_cp_h2o = '%s %s %s/%s_h2o.rda' % (cmd_handle_original, fname_rda_h2o , dir_series, name_h2o)
        run_cmd(cmd_cp_h2o, options)    
 
    # Create Control File
    if options.verbose:
        print '%%% Creating control and RAW files'
    create_control(dir_series, title_LCM, params, options)
    
    # Create raw file - MET
    cmd_bin2raw_met = '%s/%s/bin2raw %s/%s_met.rda %s/ met' % \
        (options.dir_lcmodel, lut_scanner_type[params['StationName']], dir_series, name_series, dir_series)
    run_cmd(cmd_bin2raw_met, options)

    # Create raw file - H2O

    if not options.no_h2o:
        cmd_bin2raw_h2o = '%s/%s/bin2raw %s/%s_h2o.rda %s/ h2o' % \
            (options.dir_lcmodel, lut_scanner_type[params['StationName']], dir_series, name_h2o, dir_series)    
        run_cmd(cmd_bin2raw_h2o, options)
    
    # Call LCMODEL
    if options.verbose:
        print '%%% Calling LCMODEL'
    cmd_call_LCM = '%s/bin/lcmodel < %s/control' % (options.dir_lcmodel, dir_series)
    run_cmd(cmd_call_LCM, options)
    
    # Print output
    if not options.printer=='':
        if options.verbose:
            print '%%% Printing ps output to ' + options.printer
        cmd_print = 'lpr -P %s %s/ps' % \
            (options.printer, dir_series)
        run_cmd(cmd_print, options)

    # Update Log file
    if not options.logfile=='':
        time_stamp = str(datetime.datetime.now()).split('.')[0]  # grap timestamp, but drop ms
        line_log = '%s - %s\n' % (time_stamp, dir_series)
        if options.verbose:
            print '%%% Updating logfile'
            print ' > ' + line_log
        if not options.debug:
            file_log = open('%s/%s/%s' %(options.dir_output, params['StationName'],options.logfile),'a+')
            file_log.write(line_log)
            file_log.close()

    if options.move_original:
        cmd_rm_orig = 'rmdir %s' % (dir_input_rda, )
        run_cmd(cmd_rm_orig, options)

