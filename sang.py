import argparse
from Bio import SeqIO
import pandas as pd
# Running model
from sklearn.linear_model import LogisticRegression
# Loading model
import pickle
# Required for listing files
from os import listdir
from os.path import isfile, join
import os

# Making the input robust to various 'boolean' inputs:
def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

# Converting the ab1 to the fasta:
def abi_to_seq(input_ab1_file):
    # Opening the abi file:
    test_record = SeqIO.read(input_ab1_file, 'abi')
    # Reading in the sequence:
    letters = test_record.annotations['abif_raw']['PBAS1']
    return letters
###COMBINE THE FOLLOWING TWO FUNCTIONS??
# Listing all ab1 files in directory
def listing_ab1_files(input_dir):
    # Getting all of the ab1 files:
    onlyfiles = [f for f in listdir(input_dir) if
                 isfile(join(input_dir, f)) if '.ab1' in f]
    # Throwing on the directory to the front of the ab1 filenames:
    outputlfiles = ['%s' % input_dir + each_file for each_file in onlyfiles]
    return outputlfiles, onlyfiles

# Listing all ab1 files in directory
def listing_temp_files(input_dir):
    # Getting all of the temp files:
    onlyfiles = [f for f in listdir(input_dir) if
                 isfile(join(input_dir, f)) if 'temp.ab1.conv.' in f]
    # Throwing on the directory to the front of the ab1 filenames:
    return onlyfiles

# Converting sequence to fasta file:
def seq_to_fa(input_name, input_seq, sequence_name = None):
    # Getting the sequence name for the fasta:
    if sequence_name == None:
        sequence_name = input_name

    # Generating the fasta file:
    final_filename = input_name.rsplit('.', 1)[0] + '.fa'
    final_file = open(final_filename, 'w')
    final_file.write('> %s\n' % sequence_name)
    final_file.write(input_seq + '\n')
    final_file.close()
    return final_file

# Converting ab1 file to prediction input:
def abi_to_df(input_seqio_record):
    # Reading in the abi files:
    input_seqio_record = SeqIO.read(input_seqio_record, 'abi')

    # Getting the list of letters and their locations:
    locations = list(input_seqio_record.annotations['abif_raw']['PLOC1'])
    letters = list(input_seqio_record.annotations['abif_raw']['PBAS1'])

    # Converting to df:
    letter_loc_df = pd.DataFrame()
    letter_loc_df['Locations'] = locations
    letter_loc_df['Letters'] = letters

    # Different df with all the waveform data:
    peak_df = pd.DataFrame()
    peak_df['g_let'] = list(input_seqio_record.annotations['abif_raw']['DATA9'])
    peak_df['a_let'] = list(input_seqio_record.annotations['abif_raw']['DATA10'])
    peak_df['t_let'] = list(input_seqio_record.annotations['abif_raw']['DATA11'])
    peak_df['c_let'] = list(input_seqio_record.annotations['abif_raw']['DATA12'])

    # Making the indeces play nicely and deleting the other column:
    peak_df['index_plus_one'] = peak_df.index + 1
    peak_df.index = peak_df['index_plus_one']
    letter_loc_df.index = letter_loc_df['Locations']
    letter_loc_df.drop('Locations', inplace=True, axis=1)

    # combining the dfs:
    combined_df = letter_loc_df.join(peak_df, how='inner')
    return combined_df

# Adding the previous and the following base to the df:
def surrounding_bases(input_df):
    previous_letter_value_df = input_df.shift(1)
    previous_letter_value_df.dropna(inplace=True)
    previous_letter_value_df.rename({'a_let':'prev_a','c_let':'prev_c','t_let':'prev_t','g_let':'prev_g'}, inplace=True, axis=1)

    following_letter_value_df = input_df.shift(-1)
    following_letter_value_df.dropna(inplace=True)
    following_letter_value_df.rename({'a_let':'next_a','c_let':'next_c','t_let':'next_t','g_let':'next_g'}, inplace=True, axis=1)

    current_previous_following_df = pd.concat([input_df, previous_letter_value_df, following_letter_value_df], axis=1, join='inner')
    return current_previous_following_df

def ab1_to_predicted_sequence(input_ab1_file, model):
    # Loading in and parsing input df:
    test_df = abi_to_df(input_ab1_file)
    test_letter_value_df = test_df[['a_let', 'c_let', 't_let', 'g_let']]
    test_full_info_df = surrounding_bases(test_letter_value_df)

    # Using model to predict sequence:
    predicted_probs_df = pd.DataFrame(model.predict(X=test_full_info_df),
                                      columns=['Prediction'])

    # Acquiring and returning sequence:
    sequence = ''.join(list(predicted_probs_df['Prediction']))
    return sequence

parser = argparse.ArgumentParser(description='Sanger analysis')
# Files or directories:
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-d', '--ab1_directory', metavar = '',
                    help='Directory containing ab1 files')
group.add_argument('-f', '--ab1_file', metavar = '',
                    help='ab1 file')
# If several inputs are given, have the option to make several different fastas:
parser.add_argument('-s', '--split', type=str2bool, nargs='?',
                        const=True, default=False, metavar = '',
                        help="Indicate whether to separate output fasta files")
# Name of output file
parser.add_argument('-o', '--fa_name', metavar = '',
                    help='Name of output .fa file')
# Check whether prediction is desired.
parser.add_argument('-p', '--predict', type=str2bool, nargs='?',
                        const=True, default=False, metavar = '',
                        help="Converting input sequence to predicted sequences")


args = parser.parse_args()

# Ensure either the file or the directory was selected:
if args.ab1_directory is not None and args.ab1_file is not None:
    parser.error("Please select either single file or directory of files.")

# In the event of a directory of ab1 files:
if args.ab1_directory is not None:
     # Get the file names:
     ab1_filename_list, onlyfilenames = listing_ab1_files(args.ab1_directory)

     # If one wants the predicted sequence:
     if args.predict:
         for idx, each_file in enumerate(ab1_filename_list):
             sequence = ab1_to_predicted_sequence(each_file)
             seq_to_fa('temp.ab1.conv.%s' % onlyfilenames[idx], sequence,
                        onlyfilenames[idx])

    # If one wants the assigned sequence:
     else:
        for idx, each_file in enumerate(ab1_filename_list):
            sequence = abi_to_seq(each_file)
            seq_to_fa('temp.ab1.conv.%s' % onlyfilenames[idx], sequence,
                        onlyfilenames[idx])

#### Make sure to put in the option for predicted sequence here!!!
else:
    sequence = abi_to_seq(args.ab1_file)
    seq_to_fa('temp.ab1.conv.%s' % args.ab1_file, sequence, args.ab1_file)

filenames = listing_temp_files('.')
# Making a separate file for each .fa file:
if args.split == True:
    for each_temp_file in list_of_temp_files:
        stripped_name = each_temp_file.strip('temp.ab1.conv.')
        os.rename(each_temp_file, stripped_name)

# Making a single file:
else:
    # In case a name has been set:
    if args.fa_name is not None:
        with open('%s' % args.fa_name, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())
    # In the event that a name has not been set:
    else:
        # Checking whether an file or a directory was inputted:
        if args.ab1_directory is not None:
            final_name = args.ab1_directory.replace('/', '').replace('.', '') + '.fa'
        else:
            final_name = args.ab1_file.rsplit('.', 1)[0] + '.fa'

        with open('%s' % final_name, 'w') as outfile:
            for fname in filenames:
                with open(fname) as infile:
                    outfile.write(infile.read())

# Removing all temp files:
for each_temp_file in filenames:
    os.remove(each_temp_file)
