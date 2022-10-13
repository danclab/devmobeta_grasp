import sys
import os.path as op
import mne
import pandas as pd
import scipy

from zapline_iter import zapline_until_gone

try:
    base_dir = sys.argv[1]
    pipeline = sys.argv[2]
    step_path = sys.argv[3]
except:
    print("incorrect arguments")
    sys.exit()

subjects=pd.read_csv(op.join(base_dir, 'data', 'participants.tsv'), sep='\t')
subject_ids = subjects['participant_id']

for subject_id in subject_ids:

    eeg_path = op.join(base_dir,'derivatives',pipeline,subject_id, step_path)
    epochs=mne.read_epochs_eeglab(op.join(eeg_path,'%s_task-tool_obs_exe_eeg_rereferenced_data.set' % subject_id))
    out_path = op.join(base_dir, 'derivatives', pipeline, subject_id)
    [data, iterations]=zapline_until_gone(epochs.get_data(), 60, epochs.info['sfreq'], viz=True, prefix=op.join(out_path,"zapline_iter"))

    scipy.io.savemat(op.join(out_path,'zapline.mat'),{'data':data,'iteration':iterations})
