import ConfigParser
import fnmatch
from steps import Steps


class Parameters:
    def __init__(self, initfile):

        self.cfg = ConfigParser.ConfigParser(
            {'aux_channel_source': '/virgoData/ffl/rds.ffl',
             'data_source': '/virgoData/ffl/raw.ffl'},
            allow_no_value=True)
        self.cfg.read(initfile)

        # get the BRMS parameters

        self.data_source = self.cfg.get('BRMS', 'data_source')
        self.max_freq = self.cfg.getfloat('BRMS', 'max_freq')
        self.n_points = self.cfg.getint('BRMS', 'n_points')
        self.overlap = self.cfg.getfloat('BRMS', 'overlap')
        self.steps = self.cfg.get('BRMS', 'steps').split('\n')
        self.channels_bands = {}
        steps = Steps()
        self.step_dict = {}
        for step_name, step in steps.step_dict.iteritems():
            if step_name in self.steps:
                self.step_dict[step_name] = {'class': step}
                for param in step.parameters:
                    self.step_dict[step_name][param] = self.cfg.get(step_name,
                                                                    param)

        # get the post_processing parameters
        self.aux_source = self.cfg.get('GENERAL', 'aux_channel_source')
        self.n_groups = self.cfg.getint('GENERAL', 'n_groups')
        self.all_aux = self.cfg.getboolean('GENERAL', 'all_aux')
        self.aux = self.cfg.get('GENERAL', 'aux_channels').split(
            '\n') if self.cfg.get('GENERAL', 'aux_channels') else None
        self.excluded = self.cfg.get('GENERAL', 'exclude').split(
            '\n') if self.cfg.get('GENERAL', 'exclude') else None
        self.nav = self.cfg.getint('GENERAL', 'averages')
        self.group_dict = {}
        self.aux_dict = {}
        for group_n in xrange(self.n_groups):
            sctn = 'GROUP' + str(group_n + 1)
            self.group_dict[sctn] = {'channel': self.cfg.get(sctn, 'channel'),
                                     'units': self.cfg.get(sctn, 'units'),
                                     'bands': self.cfg.get(sctn, 'bands'),
                                     'aux': self.cfg.get(sctn,
                                                         'aux_channels').split(
                                         '\n') if self.cfg.get(sctn,
                                                               'aux_channels') else None,
                                     'excl': self.cfg.get(sctn,
                                                          'exclusions').split(
                                         '\n') if self.cfg.get(sctn,
                                                               'exclusions') else None
                                     }
        for _, g_dict in self.group_dict.iteritems():
            self.extractbands(g_dict)
        if self.cfg.has_section('RESULTS'):
            self.res_param = {}
            for option in self.cfg.options('RESULTS'):

                self.res_param[option] = self.cfg.get('RESULTS', option)
                # Try to convert them to integers
                # todo: is there a better way to check if they can be converted?
                try:
                    self.res_param[option] = int(self.res_param[option])
                except ValueError:
                    pass

    # Find all the bands of the same channel, useful for the brms computation
    def merge_bands(self):
        for group, dic in self.group_dict.iteritems():
            if dic['channel'] in self.channels_bands:
                for band in dic['band_list']:
                    self.channels_bands[dic['channel']].append(band)
            else:
                self.channels_bands[dic['channel']] = dic['band_list']

    def save_extended_config(self, **kwargs):
        # todo: other interesting parameters, such as averages, npoints,
        if not self.cfg.has_section('RESULTS'):
            self.cfg.add_section('RESULTS')
        for key, arg in kwargs.iteritems():
            self.cfg.set('RESULTS', key, arg)
        with open(kwargs['hdir'] + 'config.ini', mode='w') as f:
            self.cfg.write(f)

    def extract_aux_channels(self, gpsb):
        if self.all_aux:
            chlist = get_channel_list(gpsb, self.aux_source)
            for aux in chlist:
                self.aux_dict[aux] = set(self.group_dict.keys())
        for excl_name in self.excluded or []:
            try:
                del (self.aux_dict[excl_name])
            except KeyError:
                for aux in self.aux_dict.iterkeys():
                    if (fnmatch.fnmatch(aux, 'V1:' + excl_name + '*') or
                       fnmatch.fnmatch(aux, excl_name + '*')):
                        self.excluded.append(aux)
        for aux in self.aux or []:
            self.aux_dict[aux] = set(self.group_dict.keys())
        for group_name, dic in self.group_dict.iteritems():
            for excl_name in dic['excl'] or []:
                try:
                    self.aux_dict[excl_name].remove(group_name)
                except KeyError:
                    for aux in self.aux_dict.iterkeys():
                        if (fnmatch.fnmatch(aux, 'V1:' + excl_name + '*') or
                           fnmatch.fnmatch(aux, excl_name + '*')):
                            try:
                                self.aux_dict[aux].remove(group_name)
                            except KeyError:
                                pass
            for aux in dic['aux'] or []:
                try:
                    self.aux_dict[aux].add(group_name)
                except KeyError:
                    self.aux_dict[aux] = {group_name}

    @staticmethod
    def extractbands(g_dict):
        g_dict['band_list'] = []
        b = g_dict['bands'].split(':')
        for i in range(len(b) - 1):
            g_dict['band_list'].append(b[i] + '_' + b[i + 1] + 'Hz')
