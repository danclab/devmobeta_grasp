import sys
import os.path as op
import os.path
import mne
import numpy as np
import pandas as pd
import scipy
from meegkit.dss import dss_line_iter

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


#continuous_data
#for subject_id in subject_ids:

eeg_path = op.join(base_dir, 'derivatives', pipeline, age, subject_id, '01_filtered_data')
if os.path.exists(op.join(eeg_path, '%s_task-grasping_eeg_filtered_data.set'% subject_id)):
    raw = mne.io.read_raw_eeglab(op.join(eeg_path, '%s_task-grasping_eeg_filtered_data.set'% subject_id))
    out_path = op.join(base_dir, 'derivatives', pipeline, age, subject_id)
    t_dat = raw.get_data()
    t_dat = np.transpose(t_dat, [1,0])
    [data_continuous, iterations_continuous] = dss_line_iter(t_dat, 50, raw.info['sfreq'],
                                                         prefix=op.join(out_path, "zapline_iter_cont_old"),show=True)
    data_continuous=np.transpose(data_continuous,[1,0])

    scipy.io.savemat(op.join(out_path, '%s_zapline_continuous_old.mat'% subject_id), {'data_continuous_old': data_continuous.T, 'iteration_continuous_old': iterations_continuous})
      #else:
          #continue

#epoched_data
# for subject_id in subject_ids:
#
#              eeg_path = op.join(base_dir, 'derivatives', pipeline, age, subject_id, '04_rereferenced_data')
#              print(eeg_path)
#              if os.path.exists(op.join(eeg_path, '%s_task-grasping_eeg_rereferenced_data_new.set' % subject_id)):
#                   epochs = mne.read_epochs_eeglab(op.join(eeg_path, '%s_task-grasping_eeg_rereferenced_data_new.set'% subject_id ))
#                   out_path = op.join(base_dir, 'derivatives', pipeline, age, subject_id)
#                   t_dat = epochs.get_data()
#                   t_dat = np.transpose(t_dat, [2, 1, 0])
#                   [data_epoched, iterations_epoched] = dss_line_iter(t_dat, 50, epochs.info['sfreq'],
#                                                                       prefix=op.join(out_path, "zapline_iter_epoch_new_old"),show=True)
#                   data_epoched=np.transpose(data_epoched,[2,1,0])
#                   scipy.io.savemat(op.join(out_path, '%s_zapline_epoched_new_old.mat'% subject_id ),
#                                     {'data_epoched_new_old': data_epoched.T, 'iteration_epoched_new_old': iterations_epoched})
#
#              else:
#                  continue
