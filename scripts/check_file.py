import os
import sys
def check_file(fp):
    try:
        blfile = open(fp, 'r')
    except IOError:
        blfile = open(fp, 'w')
def check_files(fp):  
    try:
        blfile = open(fp, 'r')
    except IOError:
        print("SNP file doesnot exist, download snp file to dbsnp directory, see readme")
        sys.exit()
def check_empty(ofile_path):
    if os.stat(ofile_path).st_size == 0:
        print('File is empty')
    else:
        print('File is not empty and will be overwritten')
        with open(ofile_path, 'r+') as fil:
            fil.truncate(0)