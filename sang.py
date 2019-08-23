import argparse
from Bio import SeqIO
import pandas as pd
from sklearn.linear_model import LogisticRegression
import pickle

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

parser = argparse.ArgumentParser(description='Sanger analysis')
parser.add_argument('-f', '--ab1_file', metavar = '', required=True,
                    help='Input ab1 file')
parser.add_argument('-fa', '--fa_name', metavar = '',
                    help='Name of fa file')
parser.add_argument('-p', '--predict', type=str2bool, nargs='?',
                        const=True, default=False, metavar = '',
                        help="Converting input sequence to predicted sequences")
args = parser.parse_args()

# If one wants the predicted input:
if args.predict is not None:
    # Reading in and preparing the ab1 data:
    test_df = abi_to_df(args.ab1_file)
    test_letter_value_df = test_df[['a_let', 'c_let', 't_let', 'g_let']]
    test_full_info_df = surrounding_bases(test_letter_value_df)
    # Reading in the model:
    lr = pickle.load(open('log_reg_default_million.sav', 'rb'))
    # Predicting the letters:
    predicted_probs_df = pd.DataFrame(lr.predict(X=test_full_info_df),
                                      columns=['Prediction'])
    # Getting the sequence for subsequent writing:
    sequence = ''.join(list(predicted_probs_df['Prediction']))
# If one simply wants the ab1 file letters:
else:
    sequence = abi_to_seq(args.ab1_file)

if args.fa_name is None:
    seq_to_fa(args.ab1_file, sequence)
else:
    seq_to_fa(args.fa_name, sequence)
