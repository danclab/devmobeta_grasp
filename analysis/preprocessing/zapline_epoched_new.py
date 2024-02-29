import os.path
import sys
import os.path as op
import pandas as pd
import scipy
import mne

from zapline_iter_epoched import zapline_until_gone

try:
    base_dir = sys.argv[1]
except:
    print("incorrect arguments")
    sys.exit()

subjects = pd.read_csv(op.join(base_dir, 'participants.tsv'), sep='\t')
subject_ids = subjects['participant_id']
age='twelve'
pipeline='NEARICA_NF'
subject_id='sub-204'

#for subject_id in subject_ids:

eeg_path = op.join(base_dir, 'derivatives', pipeline, age, subject_id, '04_rereferenced_data')
if os.path.exists(op.join(eeg_path, '%s_task-grasping_eeg_rereferenced_data_old.set' %(subject_id))):
    epochs = mne.read_epochs_eeglab(op.join(eeg_path, '%s_task-grasping_eeg_rereferenced_data_old.set' %(subject_id)))
    out_path = op.join(base_dir, 'derivatives', pipeline, age, subject_id)
    [data_epoched, iterations_epoched] = zapline_until_gone(epochs.get_data(), 50, epochs.info['sfreq'], viz=True,
                                                          prefix=op.join(out_path, "zapline_iter_epoch_old_new"))

    scipy.io.savemat(op.join(out_path, '%s_zapline_epoched_old_new.mat' %(subject_id)),
                                   {'data_epoched_old_new': data_epoched, 'iteration_epoched_old_new': iterations_epoched})
      #else:
          #continue
