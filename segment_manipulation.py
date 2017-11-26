import virgotools as vrg
import numpy as np
import json

# Create segments form the state vector
def check_state_vec(state_channel, gpsb, gpse, threshold,
                    condition='greater_equal'):
    if condition == 'greater_equal':
        cond_function = np.greater_equal
    else:
        if condition == 'equal':
            cond_function = np.equal
        else:
            if condition == 'greater':
                cond_function = np.greater
            else:
                if condition == 'less_equal':
                    cond_function = np.less_equal
                else:
                    if condition == 'less':
                        cond_function = np.less
                    else:
                        raise ValueError('Condition must be either '
                                         '"greater_equal", "greater",'
                                         '"equal", "less_equal" or "less"')
    print "Acquiring state vector"
    with vrg.getChannel(state_channel, gpsb, gpse - gpsb) as sttvct:
        print "Evaluating state vector and generating state flag"
        sttflag = (cond_function(sttvct.data, threshold)).astype('int')
        # finding the points where the flag changes state
        u = np.diff(np.concatenate(([0], sttflag, [0])))
        ups = np.where(u == 1, u, 0)
        # shifting the ones (raising state) by one step,
        #  removing length 1 segments
        u = (u - ups + np.roll(ups, 1))[1:]
        starts = np.argwhere(u == 1)
        stops = np.argwhere(u == -1)
        segments = np.hstack((starts, stops))/sttvct.fsample + gpsb
        print len(segments)
    return segments


def segment_subtractor(science_segments, dq_flags):
    j = 0
    for segment in science_segments:
        for dq_flag in dq_flags:
            for flagseg in dq_flag:
                if flagseg[0] < segment[0]:
                    if flagseg[1] > segment[0]:
                        if flagseg[1] < segment[1]:
                            # begin the segment after the flag ends
                            segment[0] = flagseg[1] + 1
                        else:
                            # make the length of the segment 0,
                            # it will then be removed
                            segment[0] = segment[1]
                else:
                    if flagseg[0] <= segment[1]:
                        olsegment_end = segment[1]
                        # end the segment when the flag starts
                        segment[1] = flagseg[0] - 1
                        # create a new segment from where the flag ends
                        if flagseg[1] < olsegment_end:
                            science_segments.append([flagseg[1]+1,
                                                     olsegment_end])
        j += 1


def segment_files_reader(science_segments_file, dq_flags_files):
    # Load the .json files
    dq_flags_dicts = []
    with open(science_segments_file) as f:
        science_segments_dict = json.load(f)
    for dq_flag_file in dq_flags_files:
        with open(dq_flag_file) as f:
            dq_flags_dicts.append(json.load(f))

    # TODO: check if flag known in our interval
    # Extract the active period of each flag
    active_flags = []
    science_segments = science_segments_dict['active']
    for dq_flag_dict in dq_flags_dicts:
        active_flags.append(dq_flag_dict['active'])

    # Subtract the active flag periods from the science segments
    segment_subtractor(science_segments, active_flags)

    # Find the real gpsb and gpse
    gpse = 0
    for s in science_segments:
        if s[1] == s[0]:
            science_segments.remove([s[0], s[1]])
        else:
            if s[1] > gpse:
                gpse = s[1]
    gpsb = science_segments[0][0]
    return gpsb, gpse, science_segments

