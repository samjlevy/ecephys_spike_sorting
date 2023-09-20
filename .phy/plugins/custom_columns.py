"""Show how to customize the columns in the cluster and similarity views."""

from phy import IPlugin, connect


class CustomColumnsPlugin(IPlugin):
    def attach_to_controller(self, controller):
        @connect
        def on_controller_ready(sender):
            controller.supervisor.columns = ['id', 'ch', 'depth',   \
            'amp', 'fr', 'n_spikes', 'isi_viol', 'num_viol', 'snr',  \
            'waveform_dur', 'halfwidth', 'repolarization_slope',
            'spatial_decay_slope', '%_spikes_missing', 
            'n_peaks', 'n_troughs', 'KSLabel']
            #controller.supervisor.columns = ['id', 'ch', 'depth',   \
            #'fr', 'n_spikes', 'isi_viol', 'num_viol', 'snr', 'duration',  \
            #'halfwidth', 'spread', 'presence_ratio', 'amp', 'amplitude_cutoff', \
            #'velocity_above', 'velocity_below', 'KSLabel']