
import psychopy.visual
import psychopy.event
import psychopy.clock
import psychopy.core
import psychopy.gui
import os
import pandas as pd
import numpy as np
from psychopy import core, event, parallel
import random
import datetime
import socket
import pyttsx3
import select

def send_trigger(port, code):
    # send trigger
    port.setData(code)
    # keep sending the trigger for a bit
    core.wait(0.003)
    # reset the port
    port.setData(0)
    
engine=pyttsx3.init()

font_height=60

# gui interface
gui = psychopy.gui.Dlg()
gui.addField("Subject ID:")
gui.addField("Visit")
gui.addField("Age")
gui.addField('EEG', initial=False, choices=[True, False])
gui.addField('Video Computer', initial=False, choices=[True, False])
gui.show()
print(gui.data)
if gui.OK:
    correct_input = True
    subj_id = gui.data[0]
    if not (subj_id.startswith('sub-') and subj_id.split('-')[1].isdigit() and len(subj_id.split('-')[1])==3):
        correct_input=False
        print('Bad subject ID')
    visit=gui.data[1]
    if not gui.data[1].isdigit():
         correct_input=False
         print('Bad visit number')
    else:
         visit = 'ses-{}'.format(gui.data[1])
    age = gui.data[2]
    if not gui.data[2].isdigit():
         correct_input = False
         print('Incorrect age')
    EEG = gui.data[3]
    video=gui.data[4]
else:
    core.quit()

if correct_input:

    if video:
        # connection video computer
        #TCP_IP = "192.168.2.140"
        TCP_IP = "100.1.1.2"
        TCP_PORT = 5005
        buffer_size = 1024

        # setting up the connection
        print("waiting for the video computer to init")
        s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        s.settimeout(999)
        s.bind((TCP_IP, TCP_PORT))
        s.listen(1)
        conn, addr = s.accept()

        print("CONNECTION ADDRESS:", addr)

        # waiting for the client to initialize and connect
        while True:
            try:
                data = conn.recv(buffer_size)
                if data.decode() == "connected":
                    print("video PC ready and connected")
                    break
            except:
                pass

    # file location
    df = pd.DataFrame()
    data_path = os.path.join("./data", subj_id, visit)
    if not os.path.exists(data_path):
        os.makedirs(data_path, exist_ok=True)
    timestamp=datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    logfile_fname = os.path.join(data_path, "{}_{}_{}_logfile.csv".format(subj_id, visit, timestamp))

    # win
    win = psychopy.visual.Window(fullscr=False, units="pix", color=[1, 1, 1], checkTiming=True,screen=0)

    port = parallel.ParallelPort(address=0xDFD8)
    port.setData(0)

    msg="Press a key to start event test"
    text = psychopy.visual.TextStim(win=win, text=msg, color=[-1, -1, -1],
                                    font='Arial', units='pix', height=font_height)
    text.draw()
    win.flip()
    engine.say(msg)
    engine.runAndWait()
    psychopy.event.waitKeys()

    # 1 = 'strt'
    send_trigger(port, 1)
    core.wait(1)
    # 2 = 'go'
    send_trigger(port, 2)
    core.wait(1)
    # 3 = 'grsp'
    send_trigger(port, 3)
    

    def quit_exp():
       df.to_csv(logfile_fname)
       if video:

            message_finish = "abort"
            conn.send(message_finish.encode())

            while True:
                ready = select.select([conn], [], [], 10)
                if ready[0]:
                    data = conn.recv(buffer_size)
                    if "dumped" in data.decode():
                        break

            vid_cmd = str(block).zfill(4) + "_" + "_convert_" + str(timestamp)
            conn.send(vid_cmd.encode())

            while True:
                ready = select.select([conn], [], [], 10)
                if ready[0]:
                    data = conn.recv(buffer_size)
                    if "converted" in data.decode():
                        converted_output = data.decode()
                        print(converted_output, "video data converted")
                        x, convert_time = converted_output.split("_")
                        convert_time = eval(convert_time)
                        break

            message_finish = "exit"
            conn.send(message_finish.encode())

       win.close()
       core.quit()

    event.globalKeys.add(key='q', func=quit_exp, name='shutdown')

    #callonflip
    def send_trigger_and_video_cmd(code, vid_cmd):
        if EEG:
            send_trigger(port, code)
        if video:
            conn.send(vid_cmd.encode())

    def run_block(block, block_type, block_conditions,df):

        text = psychopy.visual.TextStim(win=win, text="Block %d: %s" % ((block+1),block_type), color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
        text.draw()
        win.flip()
        engine.say('Bloc %s' % block_type)
        engine.runAndWait()
        psychopy.clock.wait(3)

        for trial in range(nTrials):
            run_trial(block, trial, block_type, block_conditions,df)
            psychopy.clock.wait(2)

        text = psychopy.visual.TextStim(win=win, text="End of Block %d.\nPress a key to continue" % ((block+1)),
                                        color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
        text.draw()
        engine.say('End of Block')
        engine.runAndWait()
        win.flip()

        engine.say('Press a key to continue')
        engine.runAndWait()

        psychopy.event.waitKeys()
        # end of bloc


    def run_trial(block, trial, block_type, block_conditions,df):

        condition = random.choice(block_conditions)


        txt="Block %d: %s\nTrial %d: %s\nPress a key when ready" % ((block+1),block_type,(trial+1),condition)
        if video and (block>0 or trial>0):
            txt="Block %d: %s\nTrial %d: %s" % ((block+1),block_type,(trial+1),condition)

        text = psychopy.visual.TextStim(win=win, text=txt, color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
        text.draw()
        win.flip()
        engine.say("%s %s" % (block_type, condition))
        engine.runAndWait()


        if (block>0 or trial>0) and video:
            dumped = False

            while True:
                if not dumped:
                    ready = select.select([conn], [], [], 10)
                    if ready[0]:
                        data = conn.recv(buffer_size)
                        if "dumped" in data.decode():
                            dump_output = data.decode()
                            print(dump_output, "video data dumped")
                            x, dump_time, x, rec_time = dump_output.split("_")
                            dump_time = float(dump_time)
                            rec_time = float(rec_time)
                            dumped=True
                            text = psychopy.visual.TextStim(win=win,
                                                            text="Block %d: %s\nTrial %d: %s.\nPress a key when ready" % ((block+1),block_type,(trial+1),condition),
                                                            color=[-1, -1, -1])
                            last_trial_row=len(df)-1
                            df.at[last_trial_row,'vid_dump_duration']=dump_time
                            df.at[last_trial_row,'vid_rec_duration']=rec_time
                            text.draw()
                            win.flip()
                            engine.say('Press a key when ready')
                            engine.runAndWait()
                keys=psychopy.event.getKeys()
                if keys is not None and len(keys) and dumped:
                    break
        else:
            text.draw()
            win.flip()
            engine.say('Press a key when ready')
            engine.runAndWait()
            while True:
                keys = psychopy.event.getKeys()
                if keys is not None and len(keys):
                    break

        video_cmd = subj_id + '_' + str(block).zfill(4) + "_" + str(trial).zfill(4) + "_start_" + str(timestamp)
        send_trigger_and_video_cmd(1, video_cmd)
        trialClock = core.Clock()
        trial_start = core.getTime()

        psychopy.clock.wait(1.5)

        #go signal
        win.callOnFlip(send_trigger, port, 2)
        text = psychopy.visual.TextStim(win=win, text="Block %d: %s\nTrial %d: %s\nShow object" % ((block+1),block_type,(trial+1),condition), color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
        text.draw()
        win.flip()
        engine.say("goooooo")
        engine.runAndWait()
        go_time = core.getTime()

        while core.getTime()-go_time<10:
            keys=psychopy.event.getKeys()
            if keys is not None and len(keys):
                break

        trial_end = trialClock.getTime()
        if keys is not None and len(keys):
            video_cmd = str(block).zfill(4) + "_" + str(trial).zfill(4) + "_stop"
            evt='grsp'
        else:
            video_cmd = "abort"
            evt = 'abrt'
        send_trigger_and_video_cmd(3, video_cmd)

        df = df.append({
            'id': subj_id,
            'visit': visit,
            'age': age,
            'block_type': block_type,
            'block': block,
            'trial': trial,
            'cond': condition,
            'trial_start': trial_start,
            'go_time': go_time,
            'trial_dur': trial_end,
            'vid_dump_duration': 0,
            'vid_rec_duration': 0
        }, ignore_index=True)
        df.to_csv(logfile_fname)


    #start
    text = psychopy.visual.TextStim(win=win, text="Start of the experiment", color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
    text.draw()
    win.flip()
    engine.say('Start of the experiment')
    engine.runAndWait()
    psychopy.clock.wait(2)
    text = psychopy.visual.TextStim(win=win, text="Press a key to continue", color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
    text.draw()
    win.flip()
    engine.say('Press a key to continue')
    engine.runAndWait()
    psychopy.event.waitKeys()

    # start of bloc
    block_types=['Cylinder','Ball']
    conditions=[['Horizontal','Vertical'],['Small','Medium','Large']]
    nBlocks = 2
    nTrials=10

    block_idx=0

    for block in range(nBlocks):

        for block_type,block_conditions in zip(block_types,conditions):
            run_block(block_idx, block_type, block_conditions,df)
            block_idx+=1

    while True:

        if event.getKeys(keyList=["c"], timeStamped=False):
            run_block(block_idx, 'Cylinder',conditions[0],df)
            block_idx += 1
        elif event.getKeys(keyList=["b"], timeStamped=False):
            run_block(block_idx, 'Ball', conditions[1],df)
            block_idx += 1


    # end of experiment
    text = psychopy.visual.TextStim(win=win, text="End of the experiment", color=[-1, -1, -1], font='Arial', units='pix', height=font_height)
    text.draw()
    engine.say('End of the experiment')
    engine.runAndWait()
    win.flip()
    psychopy.clock.wait(3, hogCPUperiod=0.2)

    quit_exp()


