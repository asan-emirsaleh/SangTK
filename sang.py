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
# Loading/running model:
import tensorflow as tf

# Running tests
import sys
import filecmp

import scipy
from scipy import stats

import numpy as np

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

def ab1_to_predicted_sequence(input_ab1_file, model, actual_ab1=True, denormalize=False):
    # Loading in and parsing input df:
    if actual_ab1 == True:
        test_df = abi_to_df(input_ab1_file)
    else:
        test_df = input_ab1_file
    test_letter_value_df = test_df[['a_let', 'c_let', 't_let', 'g_let']]
    # Rerun the nucleotide model on normalized values, but until then:
    if denormalize == True:
        test_letter_value_df = test_letter_value_df * 1000
    test_full_info_df = surrounding_bases(test_letter_value_df)

    # Using model to predict sequence:
    predicted_probs_df = pd.DataFrame(model.predict(X=test_full_info_df),
                                      columns=['Prediction'])

    # Acquiring and returning sequence:
    sequence = ''.join(list(predicted_probs_df['Prediction']))
    return sequence

# This combines a peak df with a full record
def peak_calling_df(input_df, input_seqio_record):
    input_df['peak_no_peak'] = [1] * input_df.shape[0]
    input_df.index = input_df.index + 1####MAYBE KEEP THIS IN? MAYBE REMOVE IT?
    first_val = input_df.index[0] - 5
    last_val = input_df.index[-1] + 5
    removed_df = input_df[['peak_no_peak']]
    # Different df with all the waveform data:
    peak_val = pd.DataFrame()
    peak_val['g_let'] = list(input_seqio_record.annotations['abif_raw']['DATA9'])
    peak_val['a_let'] = list(input_seqio_record.annotations['abif_raw']['DATA10'])
    peak_val['t_let'] = list(input_seqio_record.annotations['abif_raw']['DATA11'])
    peak_val['c_let'] = list(input_seqio_record.annotations['abif_raw']['DATA12'])

    peak_val = peak_val.loc[first_val:last_val]
    fin_df = removed_df.merge(peak_val, how='outer', left_index=True, right_index=True)
    zero = fin_df[fin_df['peak_no_peak'] !=1]
    zero['peak_no_peak'] = [0] * zero.shape[0]
    nonzero = fin_df[fin_df['peak_no_peak'] ==1]
    fin_df = zero.append(nonzero)
    fin_df.sort_index(inplace=True)
    return fin_df

def slope(inp_df):
    only_letters = inp_df[['g_let', 'a_let', 't_let', 'c_let']]
    slope_before = only_letters.diff(1, axis=0)
    slope_before.columns = ['slope_g_after', 'slope_a_after', 'slope_t_after', 'slope_c_after']
    slope_after = only_letters.diff(-1, axis=0)
    slope_after.columns = ['slope_g_before', 'slope_a_before', 'slope_t_before', 'slope_c_before']

    final = only_letters.join(slope_before)
    final = final.join(slope_after)
    final = final.join(inp_df[['peak_no_peak']])
    return final

def normalizing(inp_df):
    all_peak_places = inp_df[inp_df['peak_no_peak'] == 1]
    all_peak_places_vals = all_peak_places[['g_let', 'a_let', 't_let', 'c_let']]
    all_max_peaks = list(all_peak_places_vals.max(axis=1))
    trimmed_mean = scipy.stats.trim_mean(all_max_peaks, proportiontocut=0.1)
    inp_df = inp_df / trimmed_mean
    inp_df['peak_no_peak'] = inp_df['peak_no_peak'] * trimmed_mean
    inp_df['peak_no_peak'] = inp_df['peak_no_peak'].astype(int)
    return inp_df

def reshaping_the_df(inp_df, first_dim, second_dim, third_dim):
    y_val_train = np.array(inp_df['peak_no_peak'])
    inp_df = inp_df.iloc[:,:-1]
    x_val_train = inp_df.values.reshape((first_dim,second_dim,third_dim))
    return x_val_train, y_val_train

# Ensures the output sequences match previously defined output sequences
def test_check(test_kind, input_file):
    # The outputs of the nucleotide and peak sequences will be different
    if test_kind == 'nucleotide':
        validation_file = ('./test_folder/nucleotide_validation.txt')
        to_be_validated = ('./%s.fa' % input_file)
        comparison = filecmp.cmp(validation_file,to_be_validated)
        if comparison == True:
            print('Test Passed')
        else:
            print('Test Failed')

    elif test_kind == 'peak':
        validation_file = ('./test_folder/nucleotide_and_peak_validation.txt')
        to_be_validated = ('./%s.fa' % input_file)
        comparison = filecmp.cmp(validation_file,to_be_validated)
        if comparison == True:
            print('Test Passed')
        else:
            print('Test Failed')

parser = argparse.ArgumentParser(description='Sanger analysis')
# Files or directories, these are mutually exclusive options:
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-d', '--ab1_directory', metavar = '',
                    help='Directory containing ab1 files')
group.add_argument('-f', '--ab1_file', metavar = '',
                    help='ab1 file')
group.add_argument('-t', '--testing', metavar = '',
                   help='Unit and integration tests')

# If several inputs are given, have the option to make several different fastas:
parser.add_argument('-s', '--split', type=str2bool, nargs='?',
                        const=True, default=False, metavar = '',
                        help="Indicate whether to separate output fasta files")
# Name of output file
parser.add_argument('-o', '--fa_name', metavar = '',
                    help='Name of output .fa file')

# Check whether nucleotide prediction is desired.
parser.add_argument('-pn', '--predict_nucleotide', type=str2bool, nargs='?',
                        const=True, default=False, metavar = '',
                        help="Converting input sequence to predicted sequences calling nucleotides")

# Check whether peak AND nucleotide prediction is desired.
parser.add_argument('-p', '--predict_peak_and_nucleotide', type=str2bool, nargs='?',
                        const=True, default=False, metavar = '',
                        help="Converting input sequence to predicted sequences calling peaks and nucleotides")

args = parser.parse_args()

if args.testing == 'nucleotide':
    output_file = 'test_nucleotide_output'
    # Running the test for the nucleotide prediction
    os.system('python sang.py -f ab1_test_file.ab1 -o %s.fa -pn' % output_file)
    # Checking if the results match previous results
    test_check('nucleotide', output_file)
    # Cleaning up
    os.system('rm %s.fa' % output_file)
    sys.exit()
elif args.testing == 'peak':
    output_file = 'test_nucleotide_and_peak_output'
    # Running the test for the peak and nucleotide prediction
    os.system('python sang.py -f ab1_test_file.ab1 -o %s.fa -p' % output_file)
    # Checking if the results match previous results
    test_check('peak', output_file)
    # Cleaning up
    os.system('rm %s.fa' % output_file)
    sys.exit()

# Ensure either the file or the directory was selected:
if args.ab1_directory is not None and args.ab1_file is not None:
    parser.error("Please select either single file or directory of files.")

# In the event of a directory of ab1 files:
if args.ab1_directory is not None:
     # Get the file names:
     ab1_filename_list, onlyfilenames = listing_ab1_files(args.ab1_directory)

     # If one wants the nucleotide predicted sequence:
     if args.predict_nucleotide == True:
         model = pickle.load(open('log_reg_default_million.sav', 'rb'))
         for idx, each_file in enumerate(ab1_filename_list):
             sequence = ab1_to_predicted_sequence(each_file, model)
             seq_to_fa('temp.ab1.conv.%s' % onlyfilenames[idx], sequence,
                        onlyfilenames[idx])

    # If one wants the peak and the nucleotide predicted sequence:
     elif args.predict_peak_and_nucleotide == True:
        peak_model = tf.keras.models.load_model('model.h5')
        nucleotide_model = pickle.load(open('log_reg_default_million.sav', 'rb'))
        # Reading in the peak_df file(s):
        for idx, each_file in enumerate(ab1_filename_list):
            # print(each_file + '\n\n\n')
            # Reading in the file:
            current_record = SeqIO.read(each_file, 'abi')
            # Prepping it for the prediction
            current_training_df = abi_to_df(each_file)
            current_training_df['saving_og'] = current_training_df['Letters']
            fin_training = peak_calling_df(current_training_df, current_record)
            fin_training = normalizing(fin_training)
            fin_training = slope(fin_training)
            x_pred_peaks, y_pred_peaks = reshaping_the_df(fin_training,
                                                          fin_training.shape[0],
                                                          1, 12)
            # Predicting
            final_predicted_y = peak_model.predict(x_pred_peaks)
            # Getting the predicted values:
            predicted_val = []
            for each_pred in final_predicted_y:
                if each_pred[1] > 0.5:
                    predicted_val.append(1)
                else:
                    predicted_val.append(0)
            # Appending them to the dataframe:
            fin_training['predicted_peaks'] = predicted_val
            just_peaks = fin_training[fin_training['predicted_peaks']==1]

            sequence = ab1_to_predicted_sequence(just_peaks, nucleotide_model,
                                                 actual_ab1=False, denormalize=True)
            # Converting the predicted peaks to sequences:
            seq_to_fa('temp.ab1.conv.%s' % onlyfilenames[idx], sequence,
                        onlyfilenames[idx])


    # If one wants the assigned sequence:
     else:
        for idx, each_file in enumerate(ab1_filename_list):
            sequence = abi_to_seq(each_file)
            seq_to_fa('temp.ab1.conv.%s' % onlyfilenames[idx], sequence,
                        onlyfilenames[idx])

# In case you just want a single file:
else:
    # If one wants the predicted sequence:
    if args.predict_nucleotide == True:
        model = pickle.load(open('log_reg_default_million.sav', 'rb'))
        sequence = ab1_to_predicted_sequence(args.ab1_file, model)
        seq_to_fa('temp.ab1.conv.%s' % args.ab1_file, sequence, args.ab1_file)

    # If one wants the predicted peaks and nucleotides:
    elif args.predict_peak_and_nucleotide == True:
        peak_model = tf.keras.models.load_model('model.h5')
        nucleotide_model = pickle.load(open('log_reg_default_million.sav', 'rb'))
        current_record = SeqIO.read(args.ab1_file, 'abi')
        current_training_df = abi_to_df(args.ab1_file)
        current_training_df['saving_og'] = current_training_df['Letters']
        fin_training = peak_calling_df(current_training_df, current_record)
        fin_training = normalizing(fin_training)
        fin_training = slope(fin_training)
        x_pred_peaks, y_pred_peaks = reshaping_the_df(fin_training,
                                                      fin_training.shape[0],
                                                      1, 12)
        # Predicting
        final_predicted_y = peak_model.predict(x_pred_peaks)
        # Getting the predicted values:
        predicted_val = []
        for idx, item in enumerate(final_predicted_y):
            if item[1] > 0.5:
                predicted_val.append(1)
            else:
                predicted_val.append(0)
        # Appending them to the dataframe:
        fin_training['predicted_peaks'] = predicted_val
        just_peaks = fin_training[fin_training['predicted_peaks']==1]

        sequence = ab1_to_predicted_sequence(just_peaks, nucleotide_model,
                                             actual_ab1=False, denormalize=True)
        # Converting the predicted peaks to sequences:
        seq_to_fa('temp.ab1.conv.%s' % args.ab1_file, sequence, args.ab1_file)


    # If one simply wants the given sequence:
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
