from brms_config import Parameters
from segment_manipulation import check_state_vec, segment_files_reader
from BRMS_computations import ChannelParameters, process_channel
import h5py

savedir = '~/brms_new_test'


# todo: add an argument parser
def __main__(config_file, segment_params=None, segmentfiles=None):
    par = Parameters(config_file)
    par.merge_bands()
    fname = savedir + '/brms.hdf5'
    # Load the segments to use from segment files or
    # generate them from a channel and a threshold
    if segment_params:
        gpsb = segment_params[0]
        gpse = segment_params[1]
        if len(segment_params) > 2 :
            state_channel = segment_params[2]
            threshold = segment_params[3]
            try:
                condition = segment_params[4]
            except IndexError:
                condition = 'greater_equal'
            segments = check_state_vec(state_channel, gpsb, gpse, threshold,
                                          condition)
        else:
            if len(segment_params) == 2:
                segments = [(gpsb, gpse)]
            else:
                # todo: Explain better.
                raise ValueError('You should provide two or more segment '
                                 'parameters')
    else:
        if segmentfiles:
            segments = segment_files_reader(segmentfiles[0], segmentfiles[1:])
            gpsb = segments[0][0]
            gpse = segments [-1][-1]
        else:
            raise ValueError('Need to Provide segment parameters'
                             ' or segment files')

    for channel, bands in par.channels_bands.iteritems():
        ch_p = ChannelParameters(channel, par)
        results = process_channel(ch_p, par.data_source, segments)

        ########################## Data saving ###############################

        with h5py.File(fname, 'a') as f:
            f.attrs.create('gps', [gpsb, gpse])
            f.attrs.create('fsample', data=float(ch_p.ds_freq) / (ch_p.n_points * ch_p.overlap))
            f.attrs.create('segments', data=segments)
            try:
                g = f.create_group(channel)
            except ValueError:  # delete channel data if already present and create anew
                del f[channel]
                g = f.create_group(channel)
            g.create_dataset('times', data=ch_p.times)
            for band, j in zip(results[1], xrange(len(results[1]))):
                g.create_dataset('{}_{}Hz'.format(bands[j], bands[j + 1]),
                                 data=results[0][:][j])

# todo: write the resulting band in the cfg?


