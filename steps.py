import numpy as np


class Steps(object):
    def __init__(self):
        self.step_dict = {}
        for step in Steps.__subclasses__():
            self.step_dict[step.__name__] = step

    def reset(self):
        pass


class MovingMedian(Steps):
    # todo: use an abstractclass to inflict pain upon the next generations?
    # This class allows the computation of a moving median filter
    # Can be applied to both spectra and already computed brms vectors
    parameters = ['n_medians']

    def __init__(self, step_dict):
        self.buffer_size = step_dict['n_medians']
        self.median_buffer = []

    def __call__(self, pipe):
        self.median_buffer.append(pipe[1])
        if len(self.median_buffer) > self.buffer_size:
            self.median_buffer.pop(0)
        return pipe[0], np.median(self.median_buffer, axis=0)

    def reset(self):
        self.median_buffer = []


class MovingAverage(Steps):
    # todo: use an abstractclass to inflict pain upon the next generations?
    # This class allows the computation of a moving average filter
    # Can be applied to both spectra and already computed brms vectors
    parameters = ['n_averages']

    def __init__(self, step_dict):
        self.buffer_size = step_dict['n_averages']
        self.average_buffer = []
        self.sum = None

    def __call__(self, pipe):
        self.average_buffer.append(pipe[1])
        if len(self.average_buffer) > self.buffer_size:
            self.sum = np.sum((self.sum, pipe[1],
                               - self.average_buffer[0]), axis=0)
            self.average_buffer.pop(0)
        else:
            if self.sum:
                self.sum = np.sum((self.sum, pipe[1]), axis=0)
            else:
                self.sum = pipe[1]

        return pipe[0], self.sum/len(self.average_buffer)

    def reset(self):
        self.average_buffer = []
        self.sum = None


class LineRemoval(Steps):
    parameters = ['n_segments', 'n_cycles', 'n_deviations']

    def __init__(self, step_dict):
        # todo: costruisco un buffer?
        self.n_seg = step_dict['n_segments']
        self.n_cycles = step_dict['n_cycles']
        self.n_devs = step_dict['n_deviations']

    def __call__(self, pipe):
        log_spectrum = np.log(pipe[1])
        # compute the number of points for each spectrum segment/band
        n = int(np.floor(len(pipe[1]) / self.n_seg))
        # and then the actual number of segments
        nseg = int(np.floor(len(pipe[1]) / n))
        # Estimate  log spectrum floor from the minimum of each segment
        sp_floor = []
        interp_sp_floor = np.array([])
        for i in xrange(nseg):
            sp_floor.append(np.min(log_spectrum[i * n:(i + 1) * n]))
            if i != 0:
                if i == 1:
                    start_segm = np.linspace(
                        sp_floor[0] + (sp_floor[0] - sp_floor[1]) / 2,
                        sp_floor[0], np.floor(n / 2) + 1)
                    interp_sp_floor = np.concatenate((interp_sp_floor,
                                                      start_segm))
                segm = np.linspace(sp_floor[i - 1], sp_floor[i], n)
                interp_sp_floor = np.concatenate((interp_sp_floor, segm))
        end_segm = np.linspace(
            sp_floor[-1], sp_floor[1] + (sp_floor[1] - sp_floor[2]) / 2,
            np.ceil(n / 2))
        interp_sp_floor = np.concatenate((interp_sp_floor, end_segm))
        # Compute the difference between the floor and the log spectrum
        difference_log_spectrum = log_spectrum - interp_sp_floor
        # Remove the peaks from the difference log spectrum
        for j in xrange(self.n_cycles):
            # todo: There could be other ways of doing this process (i.e. using medians)
            md = np.mean(difference_log_spectrum)
            dev_thresh = self.n_devs * np.std(difference_log_spectrum)
            for i in xrange(len(difference_log_spectrum)):
                if (difference_log_spectrum[i] - md) > dev_thresh:
                    difference_log_spectrum[i] = md + dev_thresh
                if (difference_log_spectrum[i] - md) < -dev_thresh:
                    difference_log_spectrum[i] = md - dev_thresh
        # Re-estimate the log spectrum by summing
        # the peak-less difference spectrum to the floor spectrum
        nft = interp_sp_floor + difference_log_spectrum
        return pipe[0], np.exp(nft)


class BrmsComputation(Steps):
    def __init__(self, step_dict):
        self.band_idxs = (np.floor(step_dict['bands'] /
                                   step_dict['f_res'])).astype('int')
        self.n_bands = len(step_dict['bands'])

    def __call__(self, pipe):
        brmss = []
        for i in xrange(self.n_bands):
            band = pipe[1][self.band_idxs[i]:self.band_idxs[i + 1]]
            # No need to square since it already is a power spectrum
            brmss.append(np.sqrt(np.mean(band, axis=0)))
        actual_bands = pipe[0][self.band_idxs]
        return actual_bands, brmss
